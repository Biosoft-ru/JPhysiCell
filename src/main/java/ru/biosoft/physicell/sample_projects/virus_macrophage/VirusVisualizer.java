package ru.biosoft.physicell.sample_projects.virus_macrophage;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.ui.AgentVisualizer;

public class VirusVisualizer extends AgentVisualizer
{
    private int nVirus;

    public VirusVisualizer(Model model)
    {
        nVirus = model.getMicroenvironment().findDensityIndex( "virus" );
    }

    @Override
    public Color findColor(Cell pCell)
    {
        //	start with flow cytometry coloring 
        //    	std::vector<std::string> output = { "magenta" , "black" , "magenta", "black" }; 

        Color output = Color.magenta;
        //            int nVirus = microenvironment.find_density_index( "virus" );

        double min_virus = pCell.custom_data.get( "min_virion_count" );
        double max_virus = pCell.custom_data.get( "burst_virion_count" );
        double denominator = max_virus - min_virus + 1e-15;

        CellDefinition pMacrophage = CellDefinition.getCellDefinition( "macrophage" );

        // dead cells 
        if( pCell.phenotype.death.dead == true )
        {
            return Color.red.darker();
            //		 output[0] = "red"; 
            //		 output[2] = "darkred"; 
            //		 return output; 
        }

        if( pCell.type != pMacrophage.type )
        {
            output = Color.blue.darker();
            //		output[0] = "blue"; 
            //		output[2] = "darkblue"; 
            //		
            double virus = pCell.phenotype.molecular.internSubstrates[nVirus];

            if( pCell.phenotype.molecular.internSubstrates[nVirus] >= min_virus )
            {
                //                return Color.yellow;
                double interp = ( virus - min_virus ) / denominator;
                if( interp > 1.0 )
                {
                    interp = 1.0;
                }
                int Red = (int)Math.floor( 255.0 * interp );
                int Green = (int)Math.floor( 255.0 * interp );
                int Blue = (int)Math.floor( 255.0 * ( 1 - interp ) );
                return new Color( Red, Green, Blue );
            }

        }

        return output;
    }
}