package ru.biosoft.physicell.ui;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;

public class AgentColorerBW implements AgentColorer
{
    @Override
    public Color[] findColors(Cell cell)
    {
        return new Color[] {Color.black, Color.white};
    }
}
