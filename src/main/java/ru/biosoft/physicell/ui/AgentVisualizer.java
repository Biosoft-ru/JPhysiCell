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
        Color[] colors = findColors( cell );
        g.setColor( colors[0] );
        g.fillOval( x - r, y - r, 2 * r, 2 * r );
        g.setColor( colors[1] );
        g.drawOval( x - r, y - r, 2 * r, 2 * r );
    }

    public void drawAgent(Cell cell, int x, int y, int r, int nr, Graphics g)
    {
        Color[] colors = findColors( cell );
        g.setColor( colors[0] );
        g.fillOval( x - r, y - r, 2 * r, 2 * r );
        if( colors.length > 1 )
            g.setColor( colors[1] );
        else
            g.setColor( Color.black );
        g.drawOval( x - r, y - r, 2 * r, 2 * r );
        if( colors.length > 2 )
        {
            g.setColor( colors[2] );
            g.fillOval( x - nr, y - nr, 2 * nr, 2 * nr );
            g.setColor( colors[3] );
            g.drawOval( x - nr, y - nr, 2 * nr, 2 * nr );
        }
    }

    public Color findBorderColor(Cell cell)
    {
        return Color.black;
    }

    public Color findColor(Cell cell)
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

    public Color[] findColors(Cell cell)
    {
        return new Color[] {Color.black, Color.white, Color.white, Color.white};
    }

}