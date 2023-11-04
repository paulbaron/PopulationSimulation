using System;
using System.Reflection;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEditor;
using UnityEngine;

[CustomEditor(typeof(MainSimulation))]
public class MainSimulationEditor : Editor
{
    private         MainSimulation m_MainSimulation;

    void OnEnable()
    {
        if (targets.Length > 0)
            m_MainSimulation = target as MainSimulation;
        else
            m_MainSimulation = null;
    }

    void OnDisable()
    {
    }

    private void DrawAgentPropertyDistribution(float[] rawValues, float drawScale)
    {
        AnimationCurve animationCurve = new AnimationCurve();
        Keyframe[] keyframes = new Keyframe[rawValues.Length];
        
        for (int i = 0; i < keyframes.Length; ++i)
        {
            float current = rawValues[i];
            float prev = i == 0 ? current : rawValues[i - 1];
            float next = i + 1 == keyframes.Length ? current : rawValues[i + 1];
            keyframes[i].time = (float)i / (rawValues.Length - 1) * drawScale;
            keyframes[i].value = current;
            keyframes[i].inTangent = ((current - prev) * rawValues.Length) / drawScale;
            keyframes[i].outTangent = ((next - current) * rawValues.Length) / drawScale;
        }
        animationCurve.keys = keyframes;
        EditorGUILayout.CurveField(animationCurve);
    }

    private void DrawAgentPropertiesDependencies(SerializedProperty agentDependencies)
    {
        GUIStyle indent1ButtonStyle = new GUIStyle(EditorStyles.miniButton);
        GUIStyle indent2ButtonStyle = new GUIStyle(EditorStyles.miniButton);
        indent1ButtonStyle.margin = new RectOffset((EditorGUI.indentLevel + 2) * 15, 0, 0, 0);
        indent2ButtonStyle.margin = new RectOffset((EditorGUI.indentLevel + 1) * 15, 0, 0, 0);

        for (int i = 0; i < agentDependencies.arraySize; ++i)
        {
            EditorGUI.indentLevel += 1;
            SerializedProperty curDependency = agentDependencies.GetArrayElementAtIndex(i);
            SerializedProperty depPropertyFunctionOf = curDependency.FindPropertyRelative("m_PropertyFunctionOf");
            SerializedProperty depDistribution = curDependency.FindPropertyRelative("m_Distribution");
            depPropertyFunctionOf.stringValue = EditorGUILayout.TextField("Dependency Name", depPropertyFunctionOf.stringValue);
            //if (i == 0)
			{
                SerializedProperty depStdDev = curDependency.FindPropertyRelative("m_StdDev");
                depStdDev.floatValue = EditorGUILayout.FloatField("Distribution Standard Deviation", depStdDev.floatValue);
            }
            depDistribution.animationCurveValue = EditorGUILayout.CurveField(depDistribution.animationCurveValue);
            if (GUILayout.Button("Remove", indent1ButtonStyle))
            {
                agentDependencies.DeleteArrayElementAtIndex(i);
            }
            EditorGUI.indentLevel -= 1;
        }
        if (GUILayout.Button("Add Dependency", indent2ButtonStyle))
        {
            int arraySize = agentDependencies.arraySize;
            agentDependencies.InsertArrayElementAtIndex(arraySize);
            SerializedProperty newDependency = agentDependencies.GetArrayElementAtIndex(arraySize);
            SerializedProperty newDependencyPropertyFunctionOf = newDependency.FindPropertyRelative("m_PropertyFunctionOf");
            SerializedProperty newDependencyStdDev = newDependency.FindPropertyRelative("m_StdDev");
            SerializedProperty newDependencyDistribution = newDependency.FindPropertyRelative("m_Distribution");
            newDependencyPropertyFunctionOf.stringValue = "Age";
            newDependencyStdDev.floatValue = 0.12f;
            newDependencyDistribution.animationCurveValue = AnimationCurve.Constant(0, 1.0f, 1.0f);
        }
    }

    public override void OnInspectorGUI()
    {
        if (m_MainSimulation != null)
		{
            Rect rect = EditorGUILayout.GetControlRect();

            EditorGUILayout.PropertyField(serializedObject.FindProperty("m_SimulationProperties"), true);
            SerializedProperty agentProperties = serializedObject.FindProperty("m_AgentProperties");
            for (int i = 0; i < agentProperties.arraySize; ++i)
            {
                EditorGUI.indentLevel += 1;
                SerializedProperty agentProperty = agentProperties.GetArrayElementAtIndex(i);
                SerializedProperty agentPropertyName = agentProperty.FindPropertyRelative("m_Name");
                agentPropertyName.stringValue = EditorGUILayout.TextField("Property Name:", agentPropertyName.stringValue);

                SerializedProperty agentPropertyDistrib = agentProperty.FindPropertyRelative("m_PopulationDistrib");

				SerializedProperty agentPropertyDrawScale = agentProperty.FindPropertyRelative("m_DrawScale");
				SerializedProperty agentPropertyInit = agentProperty.FindPropertyRelative("m_Initialization");
				SerializedProperty agentPropertyValue = agentProperty.FindPropertyRelative("m_InitializationValue");
				SerializedProperty agentPropertyStdDev = agentProperty.FindPropertyRelative("m_InitializationStdDev");
				EditorGUILayout.PropertyField(agentPropertyDrawScale);
				EditorGUILayout.PropertyField(agentPropertyInit);
				if (agentPropertyInit.enumValueIndex == ((int)MainSimulation.InitializationFunctions.ConstantValue))
                {
                    EditorGUILayout.PropertyField(agentPropertyValue);
                }
                else if (agentPropertyInit.enumValueIndex == ((int)MainSimulation.InitializationFunctions.NormalDistribution))
                {
                    EditorGUILayout.PropertyField(agentPropertyValue);
                    EditorGUILayout.PropertyField(agentPropertyStdDev);
                }

                if (agentPropertyDistrib.arraySize != MainSimulation.BucketSize)
				{
                    agentPropertyDistrib.arraySize = MainSimulation.BucketSize;
                    for (int j = 0; j < MainSimulation.BucketSize; ++j)
                        agentPropertyDistrib.GetArrayElementAtIndex(j).floatValue = 0.0f;
                }
                float[] rawValue = new float[MainSimulation.BucketSize];
                for (int j = 0; j < MainSimulation.BucketSize; ++j)
                    rawValue[j] = agentPropertyDistrib.GetArrayElementAtIndex(j).floatValue;
                DrawAgentPropertyDistribution(rawValue, agentPropertyDrawScale.floatValue);

                SerializedProperty agentPropertyDependencies = agentProperty.FindPropertyRelative("m_Dependencies");
                DrawAgentPropertiesDependencies(agentPropertyDependencies);

                if (GUILayout.Button("Remove Property"))
                {
                    agentProperties.DeleteArrayElementAtIndex(i);
                }
                EditorGUI.indentLevel -= 1;
            }
            if (GUILayout.Button("Add Property"))
            {
                int     arraySize = agentProperties.arraySize;
                agentProperties.InsertArrayElementAtIndex(arraySize);
                SerializedProperty  newProperty = agentProperties.GetArrayElementAtIndex(arraySize);
                SerializedProperty newPropertyName = newProperty.FindPropertyRelative("m_Name");
                SerializedProperty newPropertyDistributions = newProperty.FindPropertyRelative("m_PopulationDistrib");
                newPropertyName.stringValue = "New Property";
            }
            serializedObject.ApplyModifiedProperties();
        }
    }
}