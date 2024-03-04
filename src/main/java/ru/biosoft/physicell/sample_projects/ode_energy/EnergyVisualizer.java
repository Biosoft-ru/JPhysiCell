package ru.biosoft.physicell.sample_projects.ode_energy;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.ui.FalseCellCytometryVisualizer;

public class EnergyVisualizer extends FalseCellCytometryVisualizer
{
    @Override
    public Color findColor(Cell pCell)
    {
        // Get Energy index
        int energy_vi = pCell.customData.findVariableIndex( "intra_energy" );

        // start with flow cytometry coloring 
        //	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 

        // color
        // proliferative cell
        if( pCell.phenotype.death.dead == false && pCell.type == 0 && pCell.customData.get( energy_vi ) > 445 )
        {
            return new Color( 125, 125, 0 );
            //		output[0] = "rgb(255,255,0)";
            //		output[2] = "rgb(125,125,0)";
        }

        // arrested cell
        if( pCell.phenotype.death.dead == false && pCell.type == 0 && pCell.customData.get( energy_vi ) <= 445 )
        {
            return new Color( 125, 0, 0 );
            //		output[0] = "rgb(255,0,0)";
            //		output[2] = "rgb(125,0,0)";
        }

        // dead cell
        if( pCell.phenotype.death.dead == true && pCell.type == 0 )
        {
            return new Color( 10, 10, 0 );
            //		output[0] = "rgb(20,20,20)";
            //		output[2] = "rgb(10,10,10)";
        }

        return super.findColor( pCell );//output; 
    }
}