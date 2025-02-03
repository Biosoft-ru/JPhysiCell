package ru.biosoft.physicell.ui;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;

public interface AgentColorer
{
    /**
     * returns array of colors for agent drawing
     * [0] -> agent border
     * [1] -> agent inside color
     * [2] -> agent core border color
     * [3] -> agent core inside color
     * 
     * If it returns array of length 2 then no core should be drawn (or core color is the same as inside)
     * [0] -> agent border
     * [1] -> agent inside color, agent core and agent inside color
     * 
     * If it returns array of length 1 then border and inside should have the same color
     * [0] -> agent border, inside color, core border and core inside coloer
     * @param cell
     * @return
     */
    public Color[] findColors(Cell cell);
}