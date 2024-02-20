package ru.biosoft.physicell.sample_projects.virus_macrophage;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.ui.AgentVisualizer;

public class VirusBarVisualizer extends AgentVisualizer
{
    private int nVirus;

    public VirusBarVisualizer(Model model)
    {
        nVirus = model.getMicroenvironment().findDensityIndex( "virus" );
    }

    @Override
    public Color findColor(Cell pCell)
    {
        Color output = Color.magenta;

        double min_virus = pCell.custom_data.get( "min_virion_count" );
        double max_virus = pCell.custom_data.get( "burst_virion_count" );
        double denominator = max_virus - min_virus + 1e-15;

        CellDefinition pMacrophage = CellDefinition.getCellDefinition( "macrophage" );

        // dead cells 
        if( pCell.phenotype.death.dead )
        {
            //            		 output[0] = "red"; 
            //            		 output[2] = "darkred"; 
            return Color.red.darker();
        }

        if( pCell.type != pMacrophage.type )
        {
            output = Color.blue.darker();
            //            		output[0] = "blue"; 
            //            		output[2] = "darkblue"; 

            double virus = pCell.phenotype.molecular.internalized_total_substrates[nVirus];

            if( pCell.phenotype.molecular.internalized_total_substrates[nVirus] >= min_virus )
            {
                double interp = ( virus - min_virus ) / denominator;
                if( interp > 1.0 )
                {
                    interp = 1.0;
                }

                if( interp > 0.75 )
                {
                    interp = 1;
                }
                if( interp > 0.5 && interp <= 0.75 )
                {
                    interp = 0.75;
                }
                if( interp > 0.25 && interp <= 0.5 )
                {
                    interp = 0.5;
                }
                if( interp > 0.01 && interp <= 0.25 )
                {
                    interp = 0.25;
                }

                int Red = (int)Math.floor( 255.0 * interp );
                int Green = (int)Math.floor( 255.0 * interp );
                int Blue = (int)Math.floor( 255.0 * ( 1 - interp ) );
                return new Color( Red, Green, Blue );
                //            			char szTempString [128];
                //            			sprintf( szTempString , "rgb(%u,%u,%u)", Red, Green, Blue );
                //            			output[0].assign( szTempString );
                //            			output[2].assign( szTempString );
            }

        }

        return output;
    }
}