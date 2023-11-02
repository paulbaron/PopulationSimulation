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

public class MainSimulation : MonoBehaviour
{
    public static int BucketSize = 128;
    public static float Epsilon = 1e-25f;
    public List<SimulationProperty> m_SimulationProperties = new List<SimulationProperty>();
    public List<Property> m_AgentProperties = new List<Property>();
    public List<Event> m_SimulationEvents = new List<Event>();
    public List<Event> m_AgentEvents = new List<Event>();
    public float[] m_PreviousAges = new float[BucketSize];
    public float[] m_NextAges = new float[BucketSize];

    [Serializable]
    public class SimulationProperty
    {
        public string m_Name;
        public float m_Value = 0.0f;
    }

    [Serializable]
    public class PropertyDistribution
    {
        public string m_PropertyFunctionOf;
        public AnimationCurve m_Distribution;
        public NativeArray<float> m_PrecomputeCurve;
        public float m_StdDev = 0.01f;

        public void Init()
		{
            m_PrecomputeCurve = new NativeArray<float>(BucketSize, Allocator.Persistent);
            for (int i = 0; i < BucketSize; ++i)
			{
                m_PrecomputeCurve[i] = Mathf.Clamp01(m_Distribution.Evaluate(_BucketIdxToValue(i)));
            }
		}
    }

    public class DependencyArray
    {
        public struct   RandomPick
        {
			public float m_PopulationInCell;
			public float m_RandomValue;
		}

		public NativeArray<float>   m_CurveOutput;
		public NativeArray<float>   m_PreviousPopulationDistrib = new NativeArray<float>(BucketSize, Allocator.Persistent);
        public NativeArray<float>   m_CurrentPopulationDistrib = new NativeArray<float>(BucketSize, Allocator.Persistent);
		public List<RandomPick>     m_RandomPicks = new List<RandomPick>();
		public List<RandomPick>     m_RandomPicksHistory = new List<RandomPick>();
		public NativeArray<float>   m_RestOutput;
		public float                m_RestSum = 0.0f;

		public void    Update()
        {
            // We accumulate the rest of the population per cell % 1
			m_RestOutput = new NativeArray<float>(BucketSize, Allocator.TempJob);
		    m_RestSum = 0.0f;
            float curveSum = 0.0f;
			for (int i = 0; i < BucketSize; ++i)
			{
                float rest = m_CurveOutput[i] - Mathf.Floor(m_CurveOutput[i]);
				m_CurrentPopulationDistrib[i] = Mathf.Floor(m_CurveOutput[i]);
                if (i == 0)
					m_RestOutput[i] = rest;
                else
					m_RestOutput[i] = m_RestOutput[i - 1] + rest;
				m_RestSum += rest;
                curveSum += m_CurveOutput[i];
			}
            // Keep the random picks in a list of index
            int randomPickCount = Mathf.FloorToInt(m_RestSum + 0.5f);
			// Update Random Picks
			if (randomPickCount < m_RandomPicks.Count)
			{
#if false
                // +1 on the rest of the largest values
				while (Mathf.FloorToInt(m_RestSum + 0.5f) != randomPickCount)
				{
					float maxPopInCell = 0.0f;
					int maxPopIdx = 0;
					for (int i = 0; i < BucketSize; ++i)
					{
						if (m_CurveOutput[i] >= maxPopInCell)
						{
							maxPopInCell = m_CurveOutput[i];
							maxPopIdx = i;
						}
					}
                    m_CurrentPopulationDistrib[maxPopIdx] -= 1;
					for (int i = maxPopIdx; i < BucketSize; ++i)
					{
						m_RestOutput[maxPopIdx] += 1.0f;
					}
                    m_RestSum += 1.0f;
				}
#else
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
                    m_RandomPicksHistory.Add(m_RandomPicks[maxPopIdx]);
					m_RandomPicks.RemoveAt(maxPopIdx);
    			}
#endif
			}
			else if (randomPickCount > m_RandomPicks.Count)
			{
				// Add random picks:
				while (m_RandomPicks.Count != randomPickCount)
				{
					RandomPick randomPick;
					if (m_RandomPicksHistory.Count != 0)
                    {
						randomPick = m_RandomPicksHistory[0];
                        m_RandomPicksHistory.RemoveAt(0);
					}
                    else
                    {
						randomPick = new RandomPick();

                        float randomValue = UnityEngine.Random.value * m_RestSum;
						randomPick.m_RandomValue = randomValue;
						float curCurveInt = 0.0f;

                        // Here we convert the picked sample value to a curve integral
						for (int j = 0; j < BucketSize; ++j)
						{
							curCurveInt += m_CurveOutput[j];
							if (m_RestOutput[j] >= randomValue)
							{
								float randValue = UnityEngine.Random.Range(curCurveInt - m_CurveOutput[j], curCurveInt);
								randomPick.m_RandomValue = randValue / curveSum;
								break;
							}
						}
						randomPick.m_PopulationInCell = 0.0f;
					}
					m_RandomPicks.Add(randomPick);
				}
			}
			for (int i = 0; i < m_RandomPicks.Count; ++i)
			{
                RandomPick randomPick = m_RandomPicks[i];
				float curCurveInt = 0.0f;

				for (int j = 0; j < BucketSize; ++j)
				{
                    curCurveInt += m_CurveOutput[j];
					if (curCurveInt >= randomPick.m_RandomValue * curveSum)
					{
						randomPick.m_PopulationInCell = m_CurrentPopulationDistrib[j];
						m_CurrentPopulationDistrib[j] += 1;
						break;
					}
				}
				m_RandomPicks[i] = randomPick;
			}
			m_RestOutput.Dispose();
		}
    }

	public enum InitializationFunctions
	{
        NormalDistribution,
        ConstantValue,
        RandomValue
	}

    public struct NormalDistribParallel : IJob
    {
        [ReadOnly]
        public NativeArray<float> m_PrecomputedNormalDistribution;
        [ReadOnly]
        public NativeArray<float> m_PrecomputedDependencyCurve;
        [ReadOnly]
        public NativeArray<float> m_DependencyPopulationDistrib;
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
            for (int i = m_StartIdx; i < m_StopIdx; ++i)
            {
                float curveValue = m_PrecomputedDependencyCurve[i];
                float populationDistrib = m_DependencyPopulationDistrib[i];
                int normalDistribIdxStart = _ValueToBucketIdx(curveValue - m_InvNormalPDFRangeToZero);
                int normalDistribIdxStop = _ValueToBucketIdx(curveValue + m_InvNormalPDFRangeToZero);
                int normalDistribStartIdx = (int)((1.0f - curveValue) * BucketSize);
                float previousValue = 0.0f;
                for (int k = normalDistribIdxStart; k <= normalDistribIdxStop; ++k)
                {
                    float cumulativeValue = (k == normalDistribIdxStop) ? 1.0f : m_PrecomputedNormalDistribution[normalDistribStartIdx + k];
                    float valueInBucket = cumulativeValue - previousValue;
                    previousValue = cumulativeValue;
                    m_OutputDistribution[k] += populationDistrib * valueInBucket;
                }
            }

        }
    }
	public struct DependencyNormalDistribParallel : IJob
	{
		[ReadOnly]
		public NativeArray<float> m_PrecomputedNormalDistribution;
		[ReadOnly]
		public NativeArray<float> m_PrecomputedDependencyCurve;
		[ReadOnly]
		public NativeArray<float> m_DependencyPopulationDistrib;
		[ReadOnly]
		public NativeArray<float> m_CurPropertyPopulationDistrib;
		[ReadOnly]
		public float m_PopulationCount;
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
						// Evaluate the dependency curve
						float curveValue = m_PrecomputedDependencyCurve[j];
						// Get the new propery value by multiplying the curve value with the current dependency value
						float multipliedValue = curveValue * dependencyValue;
						int normalDistribIdxStart = _ValueToBucketIdx(multipliedValue - m_InvNormalPDFRangeToZero);
						int normalDistribIdxStop = _ValueToBucketIdx(multipliedValue + m_InvNormalPDFRangeToZero);
						int normalDistribStartIdx = (int)((1.0f - multipliedValue) * BucketSize);
						// Then we redistribute the population:
						float previousValue = 0.0f;
						for (int k = normalDistribIdxStart; k <= normalDistribIdxStop; ++k)
						{
							float cumulativeValue = (k == normalDistribIdxStop) ? 1.0f : m_PrecomputedNormalDistribution[normalDistribStartIdx + k];
							float valueInBucket = cumulativeValue - previousValue;
							previousValue = cumulativeValue;
							m_OutputDistribution[k] += populationDistrib * valueInBucket;
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
        public DependencyArray m_PopulationDistrib2;
        public NativeArray<float> m_PopulationDistribNative;
		public List<PropertyDistribution> m_Dependencies = new List<PropertyDistribution>();
        public InitializationFunctions m_Initialization = InitializationFunctions.RandomValue;
        public float m_InitializationValue = 0.0f;
        public float m_InitializationStdDev = 1.0f;

        private NativeArray<float>          m_NewPopulationDistrib;
        private NativeArray<float>          m_ParallelDistrib;
		private List<NativeArray<float>>    m_ToDispose;

        public void     Destroy()
        {
            m_NewPopulationDistrib.Dispose();
            m_ParallelDistrib.Dispose();
			m_PopulationDistribNative.Dispose();
			foreach (PropertyDistribution distrib in m_Dependencies)
				distrib.m_PrecomputeCurve.Dispose();
		}

		public void     Initialize(float populationCount)
		{
			m_PopulationDistrib2 = new DependencyArray();
			m_PopulationDistrib = new float[BucketSize];
			m_ToDispose = new List<NativeArray<float>>();
			m_PopulationDistribNative = new NativeArray<float>(BucketSize, Allocator.Persistent);
			m_NewPopulationDistrib = new NativeArray<float>(BucketSize, Allocator.Persistent);
			m_ParallelDistrib = new NativeArray<float>(BucketSize * JobCount, Allocator.Persistent);

			if (m_Initialization == InitializationFunctions.RandomValue)
            {
                for (int i = 0; i < BucketSize; ++i)
					m_PopulationDistribNative[i] = populationCount / BucketSize;
            }
            else if (m_Initialization == InitializationFunctions.ConstantValue)
            {
                for (int i = 0; i < BucketSize; ++i)
					m_PopulationDistribNative[i] = 0.0f;
                float initializationValue = Mathf.Clamp01(m_InitializationValue);
				m_PopulationDistribNative[(int)(initializationValue * (BucketSize - 1))] = populationCount;
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
					m_PopulationDistribNative[i] = populationCount * valueInBucket;
                }
            }
			m_PopulationDistrib2.m_CurveOutput = m_PopulationDistribNative;
			m_PopulationDistrib2.Update();
			m_PopulationDistrib2.m_CurrentPopulationDistrib.CopyTo(m_PopulationDistrib);
			m_PopulationDistrib2.m_CurrentPopulationDistrib.CopyTo(m_PopulationDistribNative);
			for (int i = 0; i < BucketSize; ++i)
				m_NewPopulationDistrib[i] = 0;
			for (int i = 0; i < BucketSize * JobCount; ++i)
				m_ParallelDistrib[i] = 0;
			foreach (PropertyDistribution distrib in m_Dependencies)
                distrib.Init();
        }

        static int JobCount = 12;
        static readonly ProfilerMarker s_SimulatePerfMarker = new ProfilerMarker(ProfilerCategory.Ai, "PostEvaluateChange");
        public JobHandle     EvaluateChange(List<Property> properties, float populationCount, JobHandle jobDep = new JobHandle())
		{
			JobHandle lastJobDependency = jobDep;

			if (m_Dependencies.Count != 0)
			{
				NativeArray<JobHandle> handles = new NativeArray<JobHandle>(JobCount, Allocator.Temp);

				// First dependency
				{
					NativeArray<float> normalDistrib = new NativeArray<float>(BucketSize * 2, Allocator.TempJob);
					PropertyDistribution firstDependency = m_Dependencies[0];
                    float stdDev = Mathf.Max(Epsilon, firstDependency.m_StdDev);

					m_ToDispose.Add(normalDistrib);
					for (int i = 0; i < BucketSize * 2; ++i)
                    {
                        float cdfSamplingCursor = _BucketIdxToValue(i) - 1.0f;
                        float distrib = _CumulativeDistributionFunction(0.0f, stdDev, cdfSamplingCursor);
                        normalDistrib[i] = distrib;
                    }
                    // Compute a normal distribution around the curve value for the first dependency
                    Property dependentProp = properties.Find(x => x.m_Name == firstDependency.m_PropertyFunctionOf);
                    if (dependentProp != null)
                    {
                        for (int i = 0; i < JobCount; ++i)
						{
                            NormalDistribParallel job = new NormalDistribParallel();
                            job.m_PrecomputedDependencyCurve = firstDependency.m_PrecomputeCurve;
                            job.m_DependencyPopulationDistrib = dependentProp.m_PopulationDistribNative;
                            job.m_PrecomputedNormalDistribution = normalDistrib;
                            job.m_InvNormalPDFRangeToZero = _InverseNormalPDFSize(Epsilon, stdDev);
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
							//mergeJobHandle.Complete();
						}
                    }
                }
                // Multiply those values with all the other dependencies values depending on the population distribution
                // NativeArray<float> workingBuffer = new NativeArray<float>(BucketSize, Allocator.TempJob);
                for (int depIdx = 1; depIdx < m_Dependencies.Count; ++depIdx)
				{
					PropertyDistribution curDependency = m_Dependencies[depIdx];
                    Property dependentProp = properties.Find(x => x.m_Name == curDependency.m_PropertyFunctionOf);
					NativeArray<float> normalDistrib = new NativeArray<float>(BucketSize * 2, Allocator.TempJob);
					m_ToDispose.Add(normalDistrib);
					float stdDev = Mathf.Max(Epsilon, curDependency.m_StdDev);
                    for (int i = 0; i < BucketSize * 2; ++i)
                        normalDistrib[i] = 0;
                    for (int i = 0; i < BucketSize * 2; ++i)
                    {
                        float cdfSamplingCursor = _BucketIdxToValue(i) - 1.0f;
                        float distrib = _CumulativeDistributionFunction(0.0f, stdDev, cdfSamplingCursor);
                        normalDistrib[i] = distrib;
                    }
                    float invNormalPDFRangeToZero = _InverseNormalPDFSize(Epsilon, stdDev);

					for (int i = 0; i < JobCount; ++i)
					{
						DependencyNormalDistribParallel job = new DependencyNormalDistribParallel();
						job.m_PrecomputedDependencyCurve = curDependency.m_PrecomputeCurve;
						job.m_DependencyPopulationDistrib = dependentProp.m_PopulationDistribNative;
						job.m_PrecomputedNormalDistribution = normalDistrib;
						job.m_PopulationCount = populationCount;
						job.m_CurPropertyPopulationDistrib = m_NewPopulationDistrib;
						job.m_InvNormalPDFRangeToZero = _InverseNormalPDFSize(Epsilon, stdDev);
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
					}
					// workingBuffer.CopyTo(newDistribution);
				}
				handles.Dispose();
			}
			return lastJobDependency;
		}
        public void PostEvaluateChange()
        {
            s_SimulatePerfMarker.Begin();
			if (m_Dependencies.Count != 0)
            {
				float   popCount = 0.0f;
				for (int i = 0; i < BucketSize; ++i)
					popCount += m_NewPopulationDistrib[i];
				// assert ...
				foreach (NativeArray<float> buff in m_ToDispose)
					buff.Dispose();
				m_ToDispose.Clear();
				m_PopulationDistribNative.CopyFrom(m_NewPopulationDistrib);
				for (int i = 0; i < BucketSize; ++i)
					m_NewPopulationDistrib[i] = 0;
				m_PopulationDistrib2.m_CurveOutput = m_PopulationDistribNative;
				m_PopulationDistrib2.Update();
				m_PopulationDistrib2.m_CurrentPopulationDistrib.CopyTo(m_PopulationDistrib);
				m_PopulationDistrib2.m_CurrentPopulationDistrib.CopyTo(m_PopulationDistribNative);
			}
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
            m_AgentProperties[i].Initialize(populationCount.m_Value);
        Property age = m_AgentProperties.Find(x => x.m_Name == "Age");
        m_PreviousAges = new float[BucketSize];
        m_NextAges = new float[BucketSize];

        for (int i = 0; i < BucketSize; ++i)
            m_PreviousAges[i] = age.m_PopulationDistrib[i];
        Array.Copy(m_PreviousAges, 0, m_NextAges, 1, BucketSize - 1);
        m_NextAges[0] = 0;
    }

    // Update is called once per frame
    void Update()
    {
        SimulationProperty populationCount = m_SimulationProperties.Find(x => x.m_Name == "PopulationCount");
        // Update time:
        SimulationProperty speed = m_SimulationProperties.Find(x => x.m_Name == "SimulationSpeed");
        float simulationStep = Time.deltaTime * speed.m_Value;
        SimulationProperty time = m_SimulationProperties.Find(x => x.m_Name == "SimulationTime");
        float prevTimeValue = time.m_Value;
        int prevIdxValue = _ValueToBucketIdx(time.m_Value, false);
        time.m_Value += simulationStep;
        int nextIdxValue = _ValueToBucketIdx(time.m_Value, false);
        float lerpRatio = _ValueToBucketIdxMod(time.m_Value);
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
            age.m_PopulationDistrib[i] = curValue;
        }
        age.m_PopulationDistribNative.CopyFrom(age.m_PopulationDistrib);
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
                    handle = m_AgentProperties[idxToUpdate[i]].EvaluateChange(m_AgentProperties, populationCount.m_Value);
                    jobs.Add(handle);
				}
				else
                {
					handle = m_AgentProperties[idxToUpdate[i]].EvaluateChange(m_AgentProperties, populationCount.m_Value, jobs[jobs.Count - 1]);
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
				m_AgentProperties[idxToUpdate[i]].PostEvaluateChange();
			}
		}
	}
}
