package ru.biosoft.physicell.micromets;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.standard.StandardModels;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;

public class cancer_phenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );
        int apoptosis_index = phenotype.death.findDeathModelIndex( "Apoptosis" );

        phenotype.motility.isMotile = false;

        // T-cell based death
        TCell_induced_apoptosis( pCell, phenotype, dt );

        // if I am dead, remove all adhesions
        if( phenotype.death.dead == true )
        {
            // detach all attached cells
            // remove_all_adhesions( pCell );
            phenotype.secretion.secretionRates[debris_index] = pCell.customData.get(  "debris_secretion_rate");
        }

        // Mechanical contribution to proliferation
//        double mechanics_factor = strain_based_proliferation( pCell );

        // proliferation rate based on mechanical aspect
//        int cycle_G0G1_index = StandardModels.flow_cytometry_separated_cycle_model.findPhaseIndex( PhysiCellConstants.G0G1_phase );
//        int cycle_S_index =  StandardModels.flow_cytometry_separated_cycle_model.findPhaseIndex( PhysiCellConstants.S_phase );
//        pCell.phenotype.cycle.data.setTransitionRate( cycle_G0G1_index, cycle_S_index , mechanics_factor
//                * pCell.parameters.pReference_live_phenotype.cycle.data.transition_rate( cycle_G0G1_index, cycle_S_index ));//TODO:check

        // if I am dead, don't bother executing this function againtyyyyy
        if( phenotype.death.dead == true )
        {
            pCell.functions.updatePhenotype = null;
        }

        return;
    }

    void TCell_induced_apoptosis(Cell pCell, Phenotype phenotype, double dt)
    {
        int apoptosis_index = phenotype.death.findDeathModelIndex( "Apoptosis" );
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );

        if( pCell.customData.get( "TCell_contact_time" ) > pCell.customData.get( "TCell_contact_death_threshold" ) )
        {
            //            #pragma omp critical
            {
                System.out.println( "\t\t\t\t" + pCell + " (of type " + pCell.typeName + ") died from T cell contact" );
            }

            // induce death
            pCell.startDeath( apoptosis_index );
            pCell.phenotype.secretion.secretionRates[debris_index] = pCell.customData.get( "debris_secretion_rate" );

            pCell.functions.updatePhenotype = null;
        }

        return;
    }

    double strain_based_proliferation(Cell pCell)
    {
        double max_pressure = 10.0; // maximum tolerated pressure of cancer cell (proliferation) [10.0 microns]
        if( pCell.state.simplePressure < max_pressure )
        {
            return Math.pow( ( max_pressure - pCell.state.simplePressure ) / max_pressure, 1.0 );
        }
        return 0.0;
    }
}
