package ru.biosoft.physicell.sample_projects.mechano;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.ui.AgentVisualizer;

public class MechanoVisualizer extends AgentVisualizer
{
    @Override
    public Color findColor(Cell cell)
    {
        int n_springs = cell.state.springAttachments.size();
        if( !cell.typeName.equals( "BM" ) )
        {
            if( n_springs == 0 )
            {
                return Color.gray;
            }
            if( n_springs == 1 )
            {
                return new Color( 75, 0, 130 ); //indigo
            }
            if( n_springs == 2 )
            {
                return Color.blue;
            }
            if( n_springs == 3 )
            {
                return Color.green;
            }
            if( n_springs == 4 )
            {
                return Color.yellow;
            }
            if( n_springs == 5 )
            {
                return Color.orange.darker();
            }
            if( n_springs == 6 )
            {
                return Color.red;
            }
            if( n_springs > 6 )
            {
                return Color.magenta;
            }
        }
        return Color.black;
    }
}