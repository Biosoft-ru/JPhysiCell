package ru.biosoft.physicell.sample_projects.worm;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.ui.AgentVisualizer;

public class WormVisualizer extends AgentVisualizer
{
    public Color findColor(Cell cell)
    {
        if( cell.state.numberAttachedCells() == 0 )
            return Color.gray;
        //              { return { "grey", "black", "grey", "grey"}; }

        if( cell.state.numberAttachedCells() == 1 )
        {
            if( cell.customData.get( "head" ) > cell.state.attachedCells.iterator().next().customData.get( "head" ) )

                return Color.red;
            //                          return { "red", "black", "red", "red"};  } 
            return Color.orange;
        }
        //      return { "orange", "black", "orange", "orange"}; 

        //
        else if( cell.state.numberAttachedCells() >= 2 )
        {
            // shaed by head protein value 
            int intensity = (int)Math.floor( 255.0 * cell.customData.get( "head" ) );
            //                      std::string strColor = std::to_string(intensity); 
            //                      std::string color = "rgb(" + strColor + "," + strColor + ",255)"; 

            //                    if( cell.state.numberAttachedCells() > 2 )
            //                        return Color.yellow;
            //                      { return { "yellow", "black" , color, color }; }
            return new Color( intensity, intensity, 255 );
            //                      return { color , "black", color, color};    
        }
        //
        return Color.yellow;
        //  return { "yellow", "black", "yellow", "yellow" }; 
    }
}