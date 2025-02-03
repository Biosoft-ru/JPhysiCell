package ru.biosoft.physicell.ui;

import java.awt.Color;
import java.awt.Graphics;

import ru.biosoft.physicell.core.Cell;

public class AgentVisualizer
{
    private AgentColorer agentColorer = new AgentColorerDefault();

    public void drawAgent(Cell cell, int x, int y, int r, Graphics g)
    {
        Color[] colors = agentColorer.findColors( cell );
        g.setColor( colors[0] );
        g.fillOval( x - r, y - r, 2 * r, 2 * r );
        g.setColor( colors[1] );
        g.drawOval( x - r, y - r, 2 * r, 2 * r );
    }

    public void drawAgent(Cell cell, int x, int y, int r, int nr, Graphics g)
    {
        Color[] colors =  agentColorer.findColors( cell );
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

    public void setAgentColorer(AgentColorer colorer)
    {
        this.agentColorer = colorer;
    }
}