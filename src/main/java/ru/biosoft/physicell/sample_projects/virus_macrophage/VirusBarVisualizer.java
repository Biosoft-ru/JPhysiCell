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
    public Color[] findColors(Cell pCell)
    {
        Color[] output = new Color[] {Color.magenta, Color.black, Color.magenta, Color.black};
        double minVirus = pCell.customData.get( "min_virion_count" );
        double maxVirus = pCell.customData.get( "burst_virion_count" );
        double denominator = maxVirus - minVirus + 1e-15;
        CellDefinition pMacrophage = pCell.getModel().getCellDefinition( "macrophage" );
        if( pCell.phenotype.death.dead )
        {
            output[0] = Color.red;
            output[2] = Color.red.darker();
        }

        if( pCell.type != pMacrophage.type )
        {
            output[0] = Color.blue;
            output[2] = Color.blue.darker();
            double virus = pCell.phenotype.molecular.internSubstrates[nVirus];

            if( pCell.phenotype.molecular.internSubstrates[nVirus] >= minVirus )
            {
                double interp = ( virus - minVirus ) / denominator;
                if( interp > 1.0 )
                    interp = 1.0;

                if( interp > 0.75 )
                    interp = 1;

                if( interp > 0.5 && interp <= 0.75 )
                    interp = 0.75;

                if( interp > 0.25 && interp <= 0.5 )
                    interp = 0.5;

                if( interp > 0.01 && interp <= 0.25 )
                    interp = 0.25;

                int Red = (int)Math.floor( 255.0 * interp );
                int Green = (int)Math.floor( 255.0 * interp );
                int Blue = (int)Math.floor( 255.0 * ( 1 - interp ) );
                Color c = new Color( Red, Green, Blue );
                output[0] = c;
                output[2] = c;
            }
        }
        return output;
    }
}