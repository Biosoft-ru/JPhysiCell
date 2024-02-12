package ru.biosoft.physicell.sample_projects.celltypes3;

import java.awt.Color;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.ui.AgentVisualizer;

public class FluorescenceAgentVisualizer extends AgentVisualizer
{
    @Override
    public Color findColor(Cell pCell)
    {
        Microenvironment microenvironment = pCell.getMicroenvironment();
        CellDefinition aCD = CellDefinition.getCellDefinition( "A" );
        CellDefinition bCD = CellDefinition.getCellDefinition( "B" );
        CellDefinition cCD = CellDefinition.getCellDefinition( "C" );
        int A_type = aCD.type;
        int B_type = bCD.type;
        int C_type = cCD.type;

        int nA = microenvironment.findDensityIndex( "signal A" );
        int nB = microenvironment.findDensityIndex( "signal B" );
        int nC = microenvironment.findDensityIndex( "signal C" );

        // start with flow cytometry coloring 
        Color result = Color.black;
        //        String[] result = {"black", "black", "black", "black"};

        double max_fluorescence = 1; // 
        double value = 0.0;

        // color live A
        if( pCell.type == A_type )
        {
            value = pCell.phenotype.secretion.secretionRates[nA] / ( 0.001 + aCD.phenotype.secretion.secretionRates[nA] );

            value *= ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * max_fluorescence;
            if( pCell.phenotype.death.dead == true )
            {
                value = ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * max_fluorescence;
            }
            result = new Color( 1.0f, 0.0f, 1.0f, (float)value );
            //            sprintf( color, "rgba(255,0,255,%f)", value );
        }

        // color live B
        if( pCell.type == B_type )
        {
            value = pCell.phenotype.secretion.secretionRates[nB] / ( 0.001 + bCD.phenotype.secretion.secretionRates[nB] );
            value *= ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * max_fluorescence;
            if( pCell.phenotype.death.dead == true )
            {
                value = ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * max_fluorescence;
            }
            result = new Color( 0.0f, 1.0f, 0.0f, (float)value );//new Color(0,255,0);
            //            result.
            //            sprintf( color, "rgba(0,255,0,%f)", value );
        }

        // color live C
        if( pCell.type == C_type )
        {
            value = pCell.phenotype.secretion.secretionRates[nC] / ( 0.001 + cCD.phenotype.secretion.secretionRates[nC] );
            value *= ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * max_fluorescence;
            if( pCell.phenotype.death.dead == true )
            {
                value = ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * max_fluorescence;
            }
            result = new Color( 0.0f, 1.0f, 1.0f, (float)value );
            //            sprintf( color, "rgba(0,255,255,%f)", value );
        }

        // Necrotic - black
        if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_swelling
                || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_lysed
                || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic )
        {
            result = new Color( 0.0f, 0.0f, 0.0f, (float)value );
            //            sprintf( color, "rgba(0,0,0,%f)", value );
        }

        //        String[] output = {color, "none", color, "none"};
        return result;
    }
}