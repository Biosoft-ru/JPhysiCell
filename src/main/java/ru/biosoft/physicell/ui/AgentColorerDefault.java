package ru.biosoft.physicell.ui;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;

public class AgentColorerDefault implements AgentColorer
{
    public Color findBorderColor(Cell cell)
    {
        return Color.black;

    }
    public Color findColor(Cell cell)
    {
        return Color.white;
    }
    
    public Color[] findColors(Cell cell)
    {
        return new Color[] {Color.black, Color.white};
    }
}
