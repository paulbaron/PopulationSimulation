using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Profiling;
using UnityEngine.Jobs;
using Unity.Jobs;
using Unity.Collections;
using Unity.Collections.LowLevel.Unsafe;
using System.Threading;
using System.Threading.Tasks;
using System.Linq;
using System.Runtime.InteropServices.WindowsRuntime;
using UnityEditor;
using UnityEngine.Serialization;

public class MainSimulation : MonoBehaviour
{
    public static int BucketSize = 256;
    public static float Epsilon = 1e-25f;
    public List<SimulationProperty> m_SimulationProperties = new List<SimulationProperty>();
    public List<Property> m_AgentProperties = new List<Property>();
    public List<Event> m_SimulationEvents = new List<Event>();
    public List<Event> m_AgentEvents = new List<Event>();
    public float[] m_PreviousAges = new float[BucketSize];
    public float[] m_NextAges = new float[BucketSize];

	public enum PropertyType
	{
		Float,
		Int
	}

	[Serializable]
    public class SimulationProperty
    {
        public string			m_Name;
		public PropertyType		m_Type;
		public float			m_ValueF = 0.0f;
		public uint				m_ValueI = 0;
	}

	[Serializable]
    public class PropertyDistribution
    {
        public string m_PropertyFunctionOf;
        public AnimationCurve m_Distribution;
        public NativeArray<float>	m_PrecomputeCurve;
		public NativeArray<float>	m_PrecomputeNormalDistribution;
		public NativeArray<float>	m_PrecomputeNormalDistributionSum;
		public float m_StdDev = 0.01f;

        public void Init()
		{
            m_PrecomputeCurve = new NativeArray<float>(BucketSize, Allocator.Persistent);
            for (int i = 0; i < BucketSize; ++i)
                m_PrecomputeCurve[i] = Mathf.Clamp01(m_Distribution.Evaluate(_BucketIdxToValue(i)));
			m_PrecomputeNormalDistribution = new NativeArray<float>(BucketSize * 2, Allocator.Persistent);
			m_PrecomputeNormalDistributionSum = new NativeArray<float>(BucketSize * 2, Allocator.Persistent);
			float stdDev = Mathf.Max(Epsilon, m_StdDev);
			float normalDistribSum = 0.0f;
			for (int i = 0; i < BucketSize * 2; ++i)
			{
				// We could use the _CumulativeDistributionFunction here but it is less precise around 1.0f
				// And it is more intuitive to scale the distribution curve instead of having large values at minimum and maximum
				float normalDistribValue = _NormalDistributionFunction(0, stdDev, _BucketIdxToValue(i) - 1.0f);
				m_PrecomputeNormalDistribution[i] = normalDistribValue;
				normalDistribSum += normalDistribValue;
			}
			float normalDistribSumTmp = 0.0f;
			for (int i = 0; i < BucketSize * 2; ++i)
			{
				float	curNormalDistrib = m_PrecomputeNormalDistribution[i] / normalDistribSum;
				m_PrecomputeNormalDistribution[i] = curNormalDistrib;
				normalDistribSumTmp += curNormalDistrib;
				m_PrecomputeNormalDistributionSum[i] = normalDistribSumTmp;
			}
		}

		public void Destroy()
		{
			m_PrecomputeCurve.Dispose();
			m_PrecomputeNormalDistribution.Dispose();
			m_PrecomputeNormalDistributionSum.Dispose();
		}
	}

	public class PopulationDiscretizer
	{
        public struct   RandomPick
        {
			public float m_PopulationInCell;
			public float m_RandomValue;
		}

		public NativeArray<float>		m_CurveOutput;
        public NativeArray<float>		m_CurrentPopulationDistrib;
		public List<RandomPick>			m_RandomPicks = new List<RandomPick>();
		public List<RandomPick>			m_RandomPicksHistory = new List<RandomPick>();

		private NativeArray<float>	m_RestOutput;

		public void		Initialize()
		{
			m_CurrentPopulationDistrib = new NativeArray<float>(BucketSize, Allocator.Persistent);
			m_RestOutput = new NativeArray<float>(BucketSize, Allocator.Persistent);
//			m_RandomPicks = new NativeArray<RandomPick>(BucketSize, Allocator.Persistent);
		}

		public void		Destroy()
		{
			m_CurrentPopulationDistrib.Dispose();
			m_RestOutput.Dispose();
//			m_RandomPicks.Dispose();
		}

		private void	_UpdateRandomPicks(float randPickCount, float restSum)
		{
			while (m_RandomPicksHistory.Count != 0 && m_RandomPicks.Count < randPickCount)
			{
				float	currentRandom = m_RandomPicksHistory[0].m_RandomValue;
				bool	inserted = false;
				for (int i = 0; i < m_RandomPicks.Count; ++i)
				{
					if (currentRandom <= m_RandomPicks[i].m_RandomValue)
					{
						m_RandomPicks.Insert(i, m_RandomPicksHistory[0]);
						m_RandomPicksHistory.RemoveAt(0);
						inserted = true;
						break;
					}
				}
				if (!inserted)
				{
					m_RandomPicks.Add(m_RandomPicksHistory[0]);
					m_RandomPicksHistory.RemoveAt(0);
				}
			}
			// Select the largest gap between 2 samples and add random picks:
			while (m_RandomPicks.Count < randPickCount)
			{
				float previousRandValue = 0.0f;
				float largestGap = 0.0f;
				float minRandValue = 0.0f;
				float maxRandValue = 1.0f;
				int insertIdx = 0;
				for (int i = 0; i < m_RandomPicks.Count; ++i)
				{
					float currentRandValue = m_RandomPicks[i].m_RandomValue;
					float gap = currentRandValue - previousRandValue;
					if (gap >= largestGap)
					{
						minRandValue = previousRandValue;
						maxRandValue = currentRandValue;
						insertIdx = i;
						largestGap = gap;
					}
					previousRandValue = m_RandomPicks[i].m_RandomValue;
				}
				{
					float gap = 1.0f - previousRandValue;
					if (gap >= largestGap)
					{
						minRandValue = previousRandValue;
						maxRandValue = 1.0f;
						insertIdx = m_RandomPicks.Count;
						largestGap = gap;
					}
				}
				RandomPick randomPick;
				randomPick = new RandomPick();
				randomPick.m_RandomValue = UnityEngine.Random.Range(minRandValue, maxRandValue);
				randomPick.m_PopulationInCell = 0.0f;
				if (insertIdx == m_RandomPicks.Count)
					m_RandomPicks.Add(randomPick);
				else
					m_RandomPicks.Insert(insertIdx, randomPick);
			}
			// Remove random picks
			// We remove the random pick where there is the most population
			while (m_RandomPicks.Count > randPickCount)
			{
				float maxPopInCell = 0.0f;
				int maxPopIdx = 0;
				for (int i = 0; i < m_RandomPicks.Count; ++i)
				{
					if (m_RandomPicks[i].m_PopulationInCell >= maxPopInCell)
					{
						maxPopInCell = m_RandomPicks[i].m_PopulationInCell;
						maxPopIdx = i;
					}
				}
				m_RandomPicksHistory.Add(m_RandomPicks[maxPopIdx]);
				m_RandomPicks.RemoveAt(maxPopIdx);
			}
#if false
				int prevFloorValue = 0;
				for (int i = 0; i < m_RandomPicks.Count && randPickCount > m_RandomPicks.Count; ++i)
				{
					float pickValue = m_RandomPicks[i].m_RandomValue * randPickCount;
					int floorValue = Mathf.FloorToInt(pickValue + Epsilon);

					if (i < floorValue)
					{
						for (int j = 0; j < )
						float randomPos = (i + UnityEngine.Random.value) / randPickCount;
						RandomPick randomPick;
						randomPick = new RandomPick();
						randomPick.m_RandomValue = randomPos;
						randomPick.m_PopulationInCell = 0.0f;
						m_RandomPicks.Insert(i, randomPick);
						++i;
					}
//					if (i > floorValue)
//					{
//						m_RandomPicks.RemoveAt(i);
//						--i;
//					}
//					else
						prevFloorValue = floorValue;
				}
				// Add the last random picks:
				{
					int idx = Math.Max(0, m_RandomPicks.Count - 1);
					while (m_RandomPicks.Count < randPickCount)
					{
						float randomPos = (idx + UnityEngine.Random.value) / randPickCount;
						RandomPick randomPick;
						randomPick = new RandomPick();
						randomPick.m_RandomValue = randomPos;
						randomPick.m_PopulationInCell = 0.0f;
						m_RandomPicks.Add(randomPick);
						++idx;
					}
				}
			}
			if (randPickCount < m_RandomPicks.Count)
			{
				// Remove random picks
				// We remove the random pick where there is the most population
				while (m_RandomPicks.Count != randPickCount)
				{
					float maxPopInCell = 0.0f;
					int maxPopIdx = 0;
					for (int i = 0; i < m_RandomPicks.Count; ++i)
					{
						if (m_RandomPicks[i].m_PopulationInCell >= maxPopInCell)
						{
							maxPopInCell = m_RandomPicks[i].m_PopulationInCell;
							maxPopIdx = i;
						}
					}
					m_RandomPicks.RemoveAt(maxPopIdx);
				}
#endif
			}

			public void    Update(uint populationCount)
        {
		    float	restSum = 0.0f;
			uint	curPopulationCount = 0;

			for (int i = 0; i < BucketSize; ++i)
			{
				float	populationInCell = m_CurveOutput[i];
				if (Mathf.Abs(populationInCell - Mathf.Round(populationInCell)) < 0.0001f) // Float precision error
					populationInCell = Mathf.Round(populationInCell);
				curPopulationCount += (uint)populationInCell;
				float	rest = Mathf.Max(0.0f, populationInCell - Mathf.Floor(populationInCell));
				m_CurrentPopulationDistrib[i] = Mathf.Floor(populationInCell);
                if (i == 0)
					m_RestOutput[i] = rest;
                else
					m_RestOutput[i] = m_RestOutput[i - 1] + rest;
				restSum += rest;
			}
            // Keep the random picks in a list of index
            uint randomPickCount = populationCount - curPopulationCount;
			Debug.Assert(populationCount >= curPopulationCount);
			if (populationCount < curPopulationCount)
				randomPickCount = 0;
			// Update Random Picks
			_UpdateRandomPicks(randomPickCount, restSum);
#if false
			if (randomPickCount < m_RandomPicks.Count)
			{
				// Remove random picks
				// We remove the random pick where there is the most population
    			while (m_RandomPicks.Count != randomPickCount)
    			{
    				float maxPopInCell = 0.0f;
    				int maxPopIdx = 0;
    				for (int i = 0; i < m_RandomPicks.Count; ++i)
    				{
    					if (m_RandomPicks[i].m_PopulationInCell >= maxPopInCell)
    					{
    						maxPopInCell = m_RandomPicks[i].m_PopulationInCell;
                            maxPopIdx = i;
    					}
    				}
					m_RandomPicks.RemoveAt(maxPopIdx);
    			}
			}
			else if (randomPickCount > m_RandomPicks.Count)
			{
				// Add random picks:
				while (m_RandomPicks.Count != randomPickCount)
				{
					_AddRandomPick(restSum);
				}
			}
#endif
			for (int i = 0; i < m_RandomPicks.Count; ++i)
			{
                RandomPick randomPick = m_RandomPicks[i];

				for (int j = 0; j < BucketSize; ++j)
				{
					if (m_RestOutput[j] >= randomPick.m_RandomValue * restSum)
					{
						randomPick.m_PopulationInCell = m_CurrentPopulationDistrib[j];
						m_CurrentPopulationDistrib[j] += 1;
						m_RandomPicks[i] = randomPick;
						break;
					}
				}
			}
		}
    }

	public enum InitializationFunctions
	{
        NormalDistribution,
        ConstantValue,
		RandomValue,
		RandomBool
	}

	public struct NormalDistribParallel : IJob
    {
		[ReadOnly]
		public NativeArray<float>	m_PrecomputedNormalDistribution;
		[ReadOnly]
		public NativeArray<float>	m_PrecomputedNormalDistributionSum;
		[ReadOnly]
		public float				m_StdDev;
		[ReadOnly]
        public NativeArray<float>	m_PrecomputedDependencyCurve;
        [ReadOnly]
        public NativeArray<float>	m_DependencyPopulationDistrib;
		[ReadOnly]
        public float				m_InvNormalPDFRangeToZero;
        [ReadOnly]
        public int					m_StartIdx;
		[ReadOnly]
		public int					m_StopIdx;

		[NativeDisableContainerSafetyRestriction]
        public NativeSlice<float> m_OutputDistribution;

        public void Execute()
        {
            for (int i = m_StartIdx; i < m_StopIdx; ++i)
            {
                float curveValue = m_PrecomputedDependencyCurve[i];
                float populationDistrib = m_DependencyPopulationDistrib[i];
				if (populationDistrib == 0.0f)
					continue;
                int normalDistribIdxStart = _ValueToBucketIdx(curveValue - m_InvNormalPDFRangeToZero);
                int normalDistribIdxStop = _ValueToBucketIdx(curveValue + m_InvNormalPDFRangeToZero);
				if (normalDistribIdxStart == normalDistribIdxStop)
				{
					m_OutputDistribution[normalDistribIdxStart] += populationDistrib;
				}
				else
				{
					int normalDistribOffsetIdx = (int)((1.0f - curveValue) * BucketSize);
					float firstCumulativeValue = m_PrecomputedNormalDistributionSum[normalDistribOffsetIdx + normalDistribIdxStart];
					m_OutputDistribution[normalDistribIdxStart] += populationDistrib * firstCumulativeValue;
					for (int k = normalDistribIdxStart + 1; k < normalDistribIdxStop; ++k)
					{
						float currentValue = m_PrecomputedNormalDistribution[normalDistribOffsetIdx + k];
						m_OutputDistribution[k] += populationDistrib * currentValue;
					}
					float lastCumulativeValue = m_PrecomputedNormalDistributionSum[normalDistribOffsetIdx + normalDistribIdxStop - 1];
					m_OutputDistribution[normalDistribIdxStop] += populationDistrib * (1.0f - lastCumulativeValue);
				}
			}

		}
    }
	public struct DependencyNormalDistribParallel : IJob
	{
		[ReadOnly]
		public NativeArray<float>	m_PrecomputedNormalDistribution;
		[ReadOnly]
		public NativeArray<float>	m_PrecomputedNormalDistributionSum;
		[ReadOnly]
		public float				m_PopulationCount;
		[ReadOnly]
		public NativeArray<float>	m_PrecomputedDependencyCurve;
		[ReadOnly]
		public NativeArray<float>	m_DependencyPopulationDistrib;
		[ReadOnly]
		public NativeArray<float>	m_CurPropertyPopulationDistrib;
		[ReadOnly]
		public float m_InvNormalPDFRangeToZero;
		[ReadOnly]
		public int m_StartIdx;
		[ReadOnly]
		public int m_StopIdx;

		[NativeDisableContainerSafetyRestriction]
		public NativeSlice<float> m_OutputDistribution;

		public void Execute()
		{
			for (int i = 0; i < BucketSize; ++i)
			{
				float dependencyValue = _BucketIdxToValue(i);
				float distribValue = m_CurPropertyPopulationDistrib[i];
				// We compute the percentage of the total population this column represents
				float populationPercent = distribValue / m_PopulationCount;

				if (populationPercent != 0.0f)
				{
					for (int j = m_StartIdx; j < m_StopIdx; ++j)
					{
						// We then apply this percent on all buckets 
						float populationDistrib = populationPercent * m_DependencyPopulationDistrib[j];
						if (populationDistrib == 0.0f)
							continue;
						// Evaluate the dependency curve
						float curveValue = m_PrecomputedDependencyCurve[j];
						// Get the new propery value by multiplying the curve value with the current dependency value
						float multipliedValue = curveValue * dependencyValue;
						int normalDistribIdxStart = _ValueToBucketIdx(multipliedValue - m_InvNormalPDFRangeToZero);
						int normalDistribIdxStop = _ValueToBucketIdx(multipliedValue + m_InvNormalPDFRangeToZero);
						if (normalDistribIdxStart == normalDistribIdxStop)
						{
							m_OutputDistribution[normalDistribIdxStart] += populationDistrib;
						}
						else
						{
							int normalDistribOffsetIdx = (int)((1.0f - curveValue) * BucketSize);
							float firstCumulativeValue = m_PrecomputedNormalDistributionSum[normalDistribOffsetIdx + normalDistribIdxStart];
							m_OutputDistribution[normalDistribIdxStart] += populationDistrib * firstCumulativeValue;
							for (int k = normalDistribIdxStart + 1; k < normalDistribIdxStop; ++k)
							{
								float currentValue = m_PrecomputedNormalDistribution[normalDistribOffsetIdx + k];
								m_OutputDistribution[k] += populationDistrib * currentValue;
							}
							float lastCumulativeValue = m_PrecomputedNormalDistributionSum[normalDistribOffsetIdx + normalDistribIdxStop - 1];
							m_OutputDistribution[normalDistribIdxStop] += populationDistrib * (1.0f - lastCumulativeValue);
						}
					}
				}
			}

		}
	}

	public struct MergeDistribs : IJob
    {
        [ReadOnly]
        public int m_JobCount;

        public NativeArray<float> m_FinalOutput;

		public NativeArray<float> m_JobOutputs;

		public void Execute()
        {
			for (int i = 0; i < BucketSize; ++i)
			{
				float finalOutputValue = 0.0f;
				for (int j = 0; j < m_JobCount; ++j)
				{
					finalOutputValue += m_JobOutputs[j * BucketSize + i];
				}
				m_FinalOutput[i] = finalOutputValue;
			}
			for (int i = 0; i < BucketSize * m_JobCount; ++i)
				m_JobOutputs[i] = 0.0f;
		}
	}

    [Serializable]
    public class Property
    {
        public string m_Name = "NewProperty";
		public float[] m_PopulationDistrib;
        public PopulationDiscretizer m_PopulationDistribDiscretizer;
        public NativeArray<float> m_PopulationDistribNative;
		public List<PropertyDistribution> m_Dependencies = new List<PropertyDistribution>();
        public InitializationFunctions m_Initialization = InitializationFunctions.RandomValue;
        public float m_InitializationValue = 0.0f;
		public float m_InitializationStdDev = 1.0f;
		public float m_DrawScale = 1.0f;

		private NativeArray<float>          m_NewPopulationDistrib;
        private NativeArray<float>          m_ParallelDistrib;

        public void     Destroy()
        {
            m_NewPopulationDistrib.Dispose();
            m_ParallelDistrib.Dispose();
			m_PopulationDistribNative.Dispose();
			foreach (PropertyDistribution distrib in m_Dependencies)
				distrib.Destroy();
			m_PopulationDistribDiscretizer.Destroy();
		}

		public void     Initialize(uint populationCount)
		{
			m_PopulationDistrib = new float[BucketSize];
			m_PopulationDistribNative = new NativeArray<float>(BucketSize, Allocator.Persistent);
			m_NewPopulationDistrib = new NativeArray<float>(BucketSize, Allocator.Persistent);
			m_ParallelDistrib = new NativeArray<float>(BucketSize * JobCount, Allocator.Persistent);
			m_PopulationDistribDiscretizer = new PopulationDiscretizer();
			m_PopulationDistribDiscretizer.Initialize();

			if (m_Initialization == InitializationFunctions.RandomBool)
			{
				for (int i = 0; i < BucketSize; ++i)
					m_PopulationDistribNative[i] = 0.0f;
				m_PopulationDistribNative[0] = (float)populationCount / 2.0f;
				m_PopulationDistribNative[BucketSize - 1] = (float)populationCount / 2.0f;
			}
			else if (m_Initialization == InitializationFunctions.RandomValue)
			{
				for (int i = 0; i < BucketSize; ++i)
					m_PopulationDistribNative[i] = (float)populationCount / BucketSize;
			}
			else if (m_Initialization == InitializationFunctions.ConstantValue)
            {
                for (int i = 0; i < BucketSize; ++i)
					m_PopulationDistribNative[i] = 0.0f;
                float initializationValue = Mathf.Clamp01(m_InitializationValue);
				m_PopulationDistribNative[(int)(initializationValue * (BucketSize - 1))] = (float)populationCount;
            }
            else
            {
                float previousValue = 0.0f;
                for (int i = 0; i < BucketSize; ++i)
                {
                    float cdfSamplingCursor = (float)i / (float)(BucketSize - 1);
                    float cumulativeValue = cdfSamplingCursor == 1.0f ?
                                            1.0f :
                                            _CumulativeDistributionFunction(m_InitializationValue, m_InitializationStdDev, cdfSamplingCursor);
                    float valueInBucket = cumulativeValue - previousValue;

                    previousValue = cumulativeValue;
					m_PopulationDistribNative[i] = (float)populationCount * valueInBucket;
                }
            }
			m_PopulationDistribDiscretizer.m_CurveOutput = m_PopulationDistribNative;
			m_PopulationDistribDiscretizer.Update(populationCount);
			m_PopulationDistribDiscretizer.m_CurrentPopulationDistrib.CopyTo(m_PopulationDistrib);
			m_PopulationDistribDiscretizer.m_CurrentPopulationDistrib.CopyTo(m_PopulationDistribNative);
			for (int i = 0; i < BucketSize; ++i)
				m_NewPopulationDistrib[i] = 0;
			for (int i = 0; i < BucketSize * JobCount; ++i)
				m_ParallelDistrib[i] = 0;
			foreach (PropertyDistribution distrib in m_Dependencies)
                distrib.Init();
        }

        static int JobCount = 12;
        static readonly ProfilerMarker s_SimulatePerfMarker = new ProfilerMarker(ProfilerCategory.Ai, "PostEvaluateChange");
        public JobHandle     EvaluateChange(List<Property> properties, uint populationCount, JobHandle jobDep = new JobHandle())
		{
			JobHandle lastJobDependency = jobDep;

			if (m_Dependencies.Count != 0)
			{
				NativeArray<JobHandle> handles = new NativeArray<JobHandle>(JobCount, Allocator.Temp);

				// First dependency
				{
					PropertyDistribution firstDependency = m_Dependencies[0];
					// Compute a normal distribution around the curve value for the first dependency
					Property dependentProp = properties.Find(x => x.m_Name == firstDependency.m_PropertyFunctionOf);
                    if (dependentProp != null)
                    {
                        for (int i = 0; i < JobCount; ++i)
						{
                            NormalDistribParallel job = new NormalDistribParallel();
                            job.m_PrecomputedDependencyCurve = firstDependency.m_PrecomputeCurve;
                            job.m_DependencyPopulationDistrib = dependentProp.m_PopulationDistribNative;
							job.m_PrecomputedNormalDistribution = firstDependency.m_PrecomputeNormalDistribution;
							job.m_PrecomputedNormalDistributionSum = firstDependency.m_PrecomputeNormalDistributionSum;
							job.m_StdDev = Mathf.Max(Epsilon, firstDependency.m_StdDev);
                            job.m_InvNormalPDFRangeToZero = _InverseNormalPDFSize(Epsilon, job.m_StdDev);
                            job.m_StartIdx = i * (BucketSize / JobCount);
                            job.m_StopIdx = i + 1 == JobCount ? BucketSize : (i + 1) * (BucketSize / JobCount);
							job.m_OutputDistribution = m_ParallelDistrib.Slice(i * BucketSize, BucketSize);
                            handles[i] = job.Schedule(lastJobDependency);
						}
						{
                            MergeDistribs job = new MergeDistribs();
                            job.m_JobCount = JobCount;
                            job.m_JobOutputs = m_ParallelDistrib;
                            job.m_FinalOutput = m_NewPopulationDistrib;
							lastJobDependency = job.Schedule(JobHandle.CombineDependencies(handles));
							// lastJobDependency.Complete();
						}
                    }
                }
                // Multiply those values with all the other dependencies values depending on the population distribution
                // NativeArray<float> workingBuffer = new NativeArray<float>(BucketSize, Allocator.TempJob);
                for (int depIdx = 1; depIdx < m_Dependencies.Count; ++depIdx)
				{
					PropertyDistribution curDependency = m_Dependencies[depIdx];
                    Property dependentProp = properties.Find(x => x.m_Name == curDependency.m_PropertyFunctionOf);
					for (int i = 0; i < JobCount; ++i)
					{
						DependencyNormalDistribParallel job = new DependencyNormalDistribParallel();
						job.m_PrecomputedDependencyCurve = curDependency.m_PrecomputeCurve;
						job.m_DependencyPopulationDistrib = dependentProp.m_PopulationDistribNative;
						job.m_PrecomputedNormalDistribution = curDependency.m_PrecomputeNormalDistribution;
						job.m_PrecomputedNormalDistributionSum = curDependency.m_PrecomputeNormalDistributionSum;
						job.m_PopulationCount = (float)populationCount;
						job.m_InvNormalPDFRangeToZero = _InverseNormalPDFSize(Epsilon, Mathf.Max(Epsilon, curDependency.m_StdDev));
						job.m_CurPropertyPopulationDistrib = m_NewPopulationDistrib;
						job.m_StartIdx = i * (BucketSize / JobCount);
						job.m_StopIdx = i + 1 == JobCount ? BucketSize : (i + 1) * (BucketSize / JobCount);
						job.m_OutputDistribution = m_ParallelDistrib.Slice(i * BucketSize, BucketSize);
						handles[i] = job.Schedule(lastJobDependency);
					}
					{
						MergeDistribs job = new MergeDistribs();
						job.m_JobCount = JobCount;
						job.m_JobOutputs = m_ParallelDistrib;
						job.m_FinalOutput = m_NewPopulationDistrib;
						lastJobDependency = job.Schedule(JobHandle.CombineDependencies(handles));
						// lastJobDependency.Complete();
					}
				}
				handles.Dispose();
			}
			return lastJobDependency;
		}
        public void PostEvaluateChange(uint populationCount)
        {
            s_SimulatePerfMarker.Begin();
			float   popCount = 0.0f;
			for (int i = 0; i < BucketSize; ++i)
				popCount += m_NewPopulationDistrib[i];
			Debug.Assert(Mathf.Abs(popCount - (float)populationCount) < 0.5f);
			m_PopulationDistribNative.CopyFrom(m_NewPopulationDistrib);
			for (int i = 0; i < BucketSize; ++i)
				m_NewPopulationDistrib[i] = 0;
			m_PopulationDistribDiscretizer.m_CurveOutput = m_PopulationDistribNative;
			m_PopulationDistribDiscretizer.Update(populationCount);
			m_PopulationDistribDiscretizer.m_CurrentPopulationDistrib.CopyTo(m_PopulationDistrib);
			m_PopulationDistribDiscretizer.m_CurrentPopulationDistrib.CopyTo(m_PopulationDistribNative);
			s_SimulatePerfMarker.End();
		}
	}

    [Serializable]
    public class Event
	{
        public string m_Name = "NewEvent";
        public List<PropertyDistribution> m_Dependencies = new List<PropertyDistribution>();
    }
    private void    _TagPropertiesToUpdate(string propertyName, int[] updatePriority, int recIdx = 0)
	{
        for (int i = 0; i < m_AgentProperties.Count; ++i)
        {
            PropertyDistribution dependency = m_AgentProperties[i].m_Dependencies.Find(x => x.m_PropertyFunctionOf == propertyName);
            if (dependency != null)
			{
                updatePriority[i] = Math.Max(updatePriority[i], recIdx + 1);
                _TagPropertiesToUpdate(m_AgentProperties[i].m_Name, updatePriority, recIdx + 1);
            }
        }
    }
    private void    _FindPropertiesToUpdateRec(int currentPropIdx, bool[] propNeedUpdate, List<Property> propsToUpdate)
    {
        if (propNeedUpdate[currentPropIdx] && !propsToUpdate.Contains(m_AgentProperties[currentPropIdx]))
        {
            propsToUpdate.Insert(0, m_AgentProperties[currentPropIdx]);
            for (int i = 0; i < m_AgentProperties[currentPropIdx].m_Dependencies.Count; ++i)
            {
                int idx = m_AgentProperties.FindIndex(x => x.m_Name == m_AgentProperties[currentPropIdx].m_Dependencies[i].m_PropertyFunctionOf);
                if (idx >= 0)
                    _FindPropertiesToUpdateRec(idx, propNeedUpdate, propsToUpdate);
            }
        }
    }

    private void    _FindPropertiesToUpdate(bool[] propNeedUpdate, List<Property> propsToUpdate)
	{
        for (int i = 0; i < m_AgentProperties.Count; ++i)
        {
            _FindPropertiesToUpdateRec(i, propNeedUpdate, propsToUpdate);
        }
    }

    public static float _InverseNormalPDFSize(float f, float stdDev)
    {
        float x = f * stdDev * Mathf.Sqrt(2 * Mathf.PI);
        if (x <= 0.0f)
            return f > 0.0f ? 0.0f : float.PositiveInfinity;
        return stdDev * Mathf.Sqrt(-2 * Mathf.Log(x));
    }
    static private float _NormalDistributionFunction(float mean, float stdDev, float x)
    {
        float exponent = -0.5f * Mathf.Pow((x - mean) / stdDev, 2);
        float coefficient = 1.0f / (stdDev * Mathf.Sqrt(2 * Mathf.PI));
        return coefficient * Mathf.Exp(exponent);
    }

    static private float _CumulativeDistributionFunction(float mean, float stdDev, float x)
    {
        float z = (x - mean) / stdDev;
        float cdf0 = 0.5f * (1.0f + _Erf(z / Mathf.Sqrt(2)));
        return cdf0;
    }

    // Error function (complementary error function approximation)
    static private float _Erf(float x)
    {
        float a1 = 0.254829592f;
        float a2 = -0.284496736f;
        float a3 = 1.421413741f;
        float a4 = -1.453152027f;
        float a5 = 1.061405429f;
        float p = 0.3275911f;

        float sign = (x >= 0) ? 1 : -1;
        x = Mathf.Abs(x);

        float t = 1.0f / (1.0f + p * x);
        float erf = 1.0f - ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * Mathf.Exp(-x * x);
        return sign * erf;
    }
    static private int _ValueToBucketIdx(float value, bool clamp = true)
    {
        if (clamp)
		{
            if (value <= 0.0f)
                return 0;
            else if (value >= 1.0f)
                return BucketSize - 1;
        }
        return (int)(value * BucketSize);
    }
    static private float _ValueToBucketIdxMod(float value)
    {
        return (value * BucketSize) - Mathf.Floor(value * BucketSize);
    }

    static private float _BucketIdxToValue(int bucketIdx)
    {
        return BucketSize == 1 ? 0.0f : ((float)bucketIdx / (BucketSize - 1));
    }

	private void OnDisable()
	{
		foreach (Property prop in m_AgentProperties)
            prop.Destroy();
	}

	// Start is called before the first frame update
	void Start()
    {
        SimulationProperty populationCount = m_SimulationProperties.Find(x => x.m_Name == "PopulationCount");
        // Initialize Agents Distribution
        for (int i = 0; i < m_AgentProperties.Count; ++i)
            m_AgentProperties[i].Initialize(populationCount.m_ValueI);
        Property age = m_AgentProperties.Find(x => x.m_Name == "Age");
        m_PreviousAges = new float[BucketSize];
        m_NextAges = new float[BucketSize];

        for (int i = 0; i < BucketSize; ++i)
            m_PreviousAges[i] = age.m_PopulationDistribNative[i];
        Array.Copy(m_PreviousAges, 0, m_NextAges, 1, BucketSize - 1);
        m_NextAges[0] = 0;
    }

    // Update is called once per frame
    void Update()
    {
        SimulationProperty populationCount = m_SimulationProperties.Find(x => x.m_Name == "PopulationCount");
        // Update time:
        SimulationProperty speed = m_SimulationProperties.Find(x => x.m_Name == "SimulationSpeed");
        float simulationStep = Time.deltaTime * speed.m_ValueF;
        SimulationProperty time = m_SimulationProperties.Find(x => x.m_Name == "SimulationTime");
        float prevTimeValue = time.m_ValueF;
        int prevIdxValue = _ValueToBucketIdx(time.m_ValueF, false);
        time.m_ValueF += simulationStep;
        int nextIdxValue = _ValueToBucketIdx(time.m_ValueF, false);
        float lerpRatio = _ValueToBucketIdxMod(time.m_ValueF);
        if (prevIdxValue < nextIdxValue)
		{
            int idxOffset = nextIdxValue - prevIdxValue;
            Array.Copy(m_PreviousAges, 0, m_PreviousAges, idxOffset, BucketSize - idxOffset);
            Array.Clear(m_PreviousAges, 0, idxOffset);
            Array.Copy(m_PreviousAges, 0, m_NextAges, 1, BucketSize - 1);
            m_NextAges[0] = 0.0f;
        }
        Property age = m_AgentProperties.Find(x => x.m_Name == "Age");
        for (int i = 0; i < BucketSize; ++i)
		{
            float curValue = Mathf.Lerp(m_PreviousAges[i], m_NextAges[i], lerpRatio);
            age.m_PopulationDistribNative[i] = curValue;
        }

		age.m_PopulationDistribDiscretizer.m_CurveOutput = age.m_PopulationDistribNative;
		age.m_PopulationDistribDiscretizer.Update(populationCount.m_ValueI);
		age.m_PopulationDistribDiscretizer.m_CurrentPopulationDistrib.CopyTo(age.m_PopulationDistrib);
		age.m_PopulationDistribDiscretizer.m_CurrentPopulationDistrib.CopyTo(age.m_PopulationDistribNative);

		// Update Agents Properties:
		int[]  updatePriority = new int[m_AgentProperties.Count];
        for (int i = 0; i < m_AgentProperties.Count; ++i)
            updatePriority[i] = 0;
        _TagPropertiesToUpdate("Age", updatePriority);
        // Update properties:
        List<int>   idxToUpdate = new List<int>(m_AgentProperties.Count);
        for (int i = 0; i < m_AgentProperties.Count; ++i)
            idxToUpdate.Add(i);
        idxToUpdate.Sort(delegate (int a, int b)
        {
            return updatePriority[a] - updatePriority[b];
        });
        List<JobHandle>   jobs = new List<JobHandle>();

        jobs.Add(new JobHandle());
		for (int i = 0; i < idxToUpdate.Count; ++i)
		{
			if (updatePriority[idxToUpdate[i]] > 0)
			{
                JobHandle handle;
                if (updatePriority[idxToUpdate[i]] == 1) // Meaning there is only the age as dependency
                {
                    handle = m_AgentProperties[idxToUpdate[i]].EvaluateChange(m_AgentProperties, populationCount.m_ValueI);
                    jobs.Add(handle);
				}
				else
                {
					handle = m_AgentProperties[idxToUpdate[i]].EvaluateChange(m_AgentProperties, populationCount.m_ValueI, jobs[jobs.Count - 1]);
					jobs[jobs.Count - 1] = handle;
				}

			}
		}
        foreach (JobHandle job in jobs)
            job.Complete();
		for (int i = 0; i < idxToUpdate.Count; ++i)
		{
			if (updatePriority[idxToUpdate[i]] > 0)
			{
				m_AgentProperties[idxToUpdate[i]].PostEvaluateChange(populationCount.m_ValueI);
			}
		}
	}
}
