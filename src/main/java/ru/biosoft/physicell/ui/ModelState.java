package ru.biosoft.physicell.ui;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

public class ModelState
{
    private int time;
    private List<AgentState> agentStates = new ArrayList<>();

    public void setTime(int time)
    {
        this.time = time;
    }

    public int getTime()
    {
        return time;
    }

    public int getSize()
    {
        return agentStates.size();
    }

    public Iterable<AgentState> getAgents()
    {
        return agentStates;
    }

    public void addAgentState(AgentState agentState)
    {
        agentStates.add( agentState );
    }

    public static ModelState fromString(String str)
    {
        ModelState result = new ModelState();
        String[] lines = str.split( "\n" );
        for( int i =1; i<lines.length; i++ )
        {
            result.addAgentState( AgentState.fromString( lines[i] ) );
        }
        return result;
    }

    public static class AgentState
    {
        private double[] position;
        private double radius;
        private double innerRadius;
        private Color[] colors;

        public AgentState(double[] position, double r, double innerR, Color[] colors)
        {
            this.position = position;
            this.radius = r;
            this.innerRadius = innerR;
            this.colors = colors;
        }

        public static AgentState fromString(String str)
        {
            String[] parts = str.split( "\t" );
            double[] position = new double[] {Double.parseDouble( parts[0] ), Double.parseDouble( parts[1] ),
                    Double.parseDouble( parts[2] )};
            double outerRadius = Double.parseDouble( parts[3] );
            double innerRadius = Double.parseDouble( parts[4] );
            Color[] colors = null;
            if( parts.length == 9 )
            {
                colors = new Color[4];
                colors[0] = decodeColor( parts[5] );
                colors[1] = decodeColor( parts[6] );
                colors[2] = decodeColor( parts[7] );
                colors[3] = decodeColor( parts[8] );
            }
            else if( parts.length == 7 )
            {
                colors = new Color[2];
                colors[0] = decodeColor( parts[5] );
                colors[1] = decodeColor( parts[6] );
            }
            else if( parts.length == 6 )
            {
                colors = new Color[1];
                colors[0] = decodeColor( parts[5] );
            }
            return new AgentState( position, outerRadius, innerRadius, colors );
        }

        public double[] getPosition()
        {
            return position;
        }

        public Color[] getColors()
        {
            return colors;
        }

        public double getRadius()
        {
            return radius;
        }

        public double getInnerRadius()
        {
            return innerRadius;
        }
    }

    public static Color decodeColor(String s)
    {
        s = s.substring( 1, s.length() - 1 );
        String[] parts = s.split( "," );
        int r = Integer.parseInt( parts[0].trim() );
        int g = Integer.parseInt( parts[1].trim() );
        int b = Integer.parseInt( parts[2].trim() );
        return new Color( r, g, b );
    }
}
