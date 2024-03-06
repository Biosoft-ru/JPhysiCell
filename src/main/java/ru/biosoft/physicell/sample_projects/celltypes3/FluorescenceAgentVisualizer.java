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
    public Color[] findColors(Cell pCell)
    {
        Microenvironment m = pCell.getMicroenvironment();
        CellDefinition aCD = CellDefinition.getCellDefinition( "A" );
        CellDefinition bCD = CellDefinition.getCellDefinition( "B" );
        CellDefinition cCD = CellDefinition.getCellDefinition( "C" );
        int typeA = aCD.type;
        int typeB = bCD.type;
        int typeC = cCD.type;

        int nA = m.findDensityIndex( "signal A" );
        int nB = m.findDensityIndex( "signal B" );
        int nC = m.findDensityIndex( "signal C" );

        //        Color[] result = new Color[] {Color.black, Color.black, Color.black, Color.black};

        double maxFluorescence = 1; // 
        double value = 0.0;
        Color c = Color.black;
        // color live A
        if( pCell.type == typeA )
        {
            value = pCell.phenotype.secretion.secretionRates[nA] / ( 0.001 + aCD.phenotype.secretion.secretionRates[nA] );
            value *= ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * maxFluorescence;
            if( pCell.phenotype.death.dead == true )
                value = ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * maxFluorescence;

            c = new Color( 1.0f, 0.0f, 1.0f, (float)value );
            //            sprintf( color, "rgba(255,0,255,%f)", value );
        }

        // color live B
        if( pCell.type == typeB )
        {
            value = pCell.phenotype.secretion.secretionRates[nB] / ( 0.001 + bCD.phenotype.secretion.secretionRates[nB] );
            value *= ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * maxFluorescence;
            if( pCell.phenotype.death.dead == true )
                value = ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * maxFluorescence;
            c = new Color( 0.0f, 1.0f, 0.0f, (float)value );//new Color(0,255,0);
            //            result.
            //            sprintf( color, "rgba(0,255,0,%f)", value );
        }

        // color live C
        if( pCell.type == typeC )
        {
            value = pCell.phenotype.secretion.secretionRates[nC] / ( 0.001 + cCD.phenotype.secretion.secretionRates[nC] );
            value *= ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * maxFluorescence;
            if( pCell.phenotype.death.dead == true )
                value = ( 1.0 - pCell.phenotype.volume.fluid_fraction ) * maxFluorescence;
            c = new Color( 0.0f, 1.0f, 1.0f, (float)value );
            //            sprintf( color, "rgba(0,255,255,%f)", value );
        }

        // Necrotic - black
        if( pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_swelling
                || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic_lysed
                || pCell.phenotype.cycle.currentPhase().code == PhysiCellConstants.necrotic )
        {
            c = new Color( 0.0f, 0.0f, 0.0f, (float)value );
            //            sprintf( color, "rgba(0,0,0,%f)", value );
        }

        //        String[] output = {color, "none", color, "none"};
        return new Color[] {c, c, c, c};//result;
    }
}