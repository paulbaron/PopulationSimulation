using UnityEngine;
using System.Collections;

public class SimulationGUI : MonoBehaviour
{
	public MainSimulation m_Simulation;
	public Texture2D	m_Texture;
	public float m_WindowWidth = 500.0f;
	public float m_WindowHeight = 400.0f;
	private Rect m_WindowRect;
	private float m_DisplayScale = 0.01f;
	private int m_PropertyIndex = 0;

	private void Start()
	{
		m_WindowRect = new Rect(Screen.width / 2 - m_WindowWidth / 2, Screen.height / 2 - m_WindowHeight / 2, m_WindowWidth, m_WindowHeight);
	}

	void OnGUI()
	{
		m_WindowRect = GUI.Window(0, m_WindowRect, DoMyWindow, "My Window");
	}

	// Make the contents of the window
	void DoMyWindow(int windowID)
	{
		float	graphXPos = 40.0f;
		float	graphYPos = 140.0f;
		float	graphWidth = 300.0f;
		float	graphHeight = 150.0f;

		float	displayScale2 = m_DisplayScale * m_DisplayScale;
		GUI.Box(new Rect(graphXPos - 20, graphYPos - 20, graphWidth + 40, graphHeight + 40), "Graph");
		if (m_PropertyIndex < m_Simulation.m_AgentProperties.Count)
		{
			float	colWidth = graphWidth / MainSimulation.BucketSize;
			for (int i = 0; i < MainSimulation.BucketSize; ++i)
			{
				float populationCount = m_Simulation.m_AgentProperties[m_PropertyIndex].m_PopulationDistrib[i];
				float height = Mathf.Min((populationCount * displayScale2) * graphHeight, graphHeight);
				float offsetY = graphHeight - height;
				Rect rect = new Rect(graphXPos + i * colWidth, graphYPos + offsetY, colWidth + 1, height);
				GUI.DrawTexture(rect, m_Texture);
			}
		}
		GUI.Label(new Rect(graphXPos, 20, 200, 50), "Display Scale: ");
		m_DisplayScale = GUI.HorizontalSlider(new Rect(graphXPos + 100, 20, 200, 50), m_DisplayScale, 0.0000001f, 1.0f);
		GUI.Label(new Rect(graphXPos, 70, 200, 50), "Property: ");
		int propertyIdx = (int)GUI.HorizontalSlider(new Rect(graphXPos + 100, 70, 200, 50), (float)m_PropertyIndex, 0.0f, (float)m_Simulation.m_AgentProperties.Count - 1);
		string propertyName = "";
		if (m_PropertyIndex < m_Simulation.m_AgentProperties.Count)
			propertyName = m_Simulation.m_AgentProperties[m_PropertyIndex].m_Name;
		GUI.Label(new Rect(graphXPos, 90, 200, 50), propertyName);
		if (m_PropertyIndex != propertyIdx)
		{
			m_PropertyIndex = propertyIdx;
			// Rescale the display size:
			if (m_PropertyIndex < m_Simulation.m_AgentProperties.Count)
			{
				float maxPopulationCount = 0.0f;
				for (int i = 0; i < MainSimulation.BucketSize; ++i)
				{
					float populationCount = m_Simulation.m_AgentProperties[m_PropertyIndex].m_PopulationDistrib[i];
					maxPopulationCount = Mathf.Max(maxPopulationCount, populationCount);
				}
				m_DisplayScale = Mathf.Sqrt(1.0f / maxPopulationCount);
			}
		}
		GUI.DragWindow(new Rect(0, 0, 10000, 10000));
	}
}