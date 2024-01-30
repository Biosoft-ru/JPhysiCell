package ru.biosoft.physicell.ui;

import java.awt.Color;
import java.awt.Graphics;
import java.util.HashMap;
import java.util.Map;

import ru.biosoft.physicell.core.Cell;

public class AgentVisualizer
{
    Map<Integer, Color> typeColor = new HashMap<>();
    Map<String, Color> phaseColor = new HashMap<>();

    public void addTypeColor(Integer type, Color color)
    {
        typeColor.put( type, color );
    }

    public void addPhaseColor(String type, Color color)
    {
        phaseColor.put( type, color );
    }

    public void drawAgent(Cell cell, int x, int y, int r, Graphics g)
    {
        g.setColor( findColor( cell ) );
        g.fillOval( x - r, y - r, 2 * r, 2 * r );
        g.setColor( Color.black );
        g.drawOval( x - r, y - r, 2 * r, 2 * r );
    }

    private Color findColor(Cell cell)
    {
        Color c = Color.white;
        int type = cell.type;
        if( typeColor.containsKey( type ) )
            c = typeColor.get( type );
        else
        {
            String phase = cell.phenotype.cycle.currentPhase().name;
            if( phaseColor.containsKey( phase ) )
                c = phaseColor.get( phase );
        }
        return c;
    }

}