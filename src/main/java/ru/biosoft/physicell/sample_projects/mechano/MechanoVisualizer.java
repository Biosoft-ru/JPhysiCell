package ru.biosoft.physicell.sample_projects.mechano;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.ui.AgentColorer;

public class MechanoVisualizer implements AgentColorer
{
    @Override
    public Color[] findColors(Cell cell)
    {
        int springs = cell.state.springAttachments.size();
        Color[] result = new Color[] {Color.black};
        if( !cell.typeName.equals( "BM" ) )
        {
            if( springs == 0 )
                result[0] = Color.gray;
            else if( springs == 1 )
                result[0] = new Color( 75, 0, 130 ); //indigo
            else if( springs == 2 )
                result[0] = Color.blue;
            else if( springs == 3 )
                result[0] = Color.green;
            else if( springs == 4 )
                result[0] = Color.yellow;
            else if( springs == 5 )
                result[0] = Color.orange.darker();
            else if( springs == 6 )
                result[0] = Color.red;
            else if( springs > 6 )
                result[0] = Color.magenta;
        }
        return result;
    }
}