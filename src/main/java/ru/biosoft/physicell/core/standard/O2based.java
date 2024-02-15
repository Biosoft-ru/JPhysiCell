package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;

public class O2based extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        // supported cycle models:
        // advanced_Ki67_cycle_model= 0;
        // basic_Ki67_cycle_model=1
        // live_cells_cycle_model = 5; 

        if( phenotype.death.dead == true )
            return;

        // set up shortcuts to find the Q and K(1) phases (assuming Ki67 basic or advanced model)
        boolean indices_initiated = false;
        int start_phase_index = 0; // Q_phase_index; 
        int end_phase_index = 0; // K_phase_index;
        int necrosis_index = 0;

        int oxygen_substrate_index = pCell.getMicroenvironment().findDensityIndex( "oxygen" );

        if( !indices_initiated )
        {
            // Ki67 models
            if( phenotype.cycle.code == PhysiCellConstants.advanced_Ki67_cycle_model
                    || phenotype.cycle.code == PhysiCellConstants.basic_Ki67_cycle_model )
            {
                start_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.Ki67_negative );
                necrosis_index = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );

                if( phenotype.cycle.code == PhysiCellConstants.basic_Ki67_cycle_model )
                {
                    end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.Ki67_positive );
                    indices_initiated = true;
                }
                if( phenotype.cycle.code == PhysiCellConstants.advanced_Ki67_cycle_model )
                {
                    end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.Ki67_positive_premitotic );
                    indices_initiated = true;
                }
            }

            // live model 
            if( phenotype.cycle.code == PhysiCellConstants.live_cells_cycle_model )
            {
                start_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.live );
                necrosis_index = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );
                end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.live );
                indices_initiated = true;
            }

            // cytometry models 
            if( phenotype.cycle.code == PhysiCellConstants.flow_cytometry_cycle_model
                    || phenotype.cycle.code == PhysiCellConstants.flow_cytometry_separated_cycle_model )
            {
                start_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.G0G1_phase );
                necrosis_index = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );
                end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.S_phase );
                indices_initiated = true;
            }

            if( phenotype.cycle.code == PhysiCellConstants.cycling_quiescent_model )
            {
                start_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.quiescent );
                necrosis_index = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );
                end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.cycling );
                indices_initiated = true;
            }

        }

        // don't continue if we never "figured out" the current cycle model. 
        if( !indices_initiated )
            return;

        // sample the microenvironment to get the pO2 value 
        double pO2 = ( pCell.nearest_density_vector() )[oxygen_substrate_index]; // PhysiCellConstants.oxygen_index]; 
        int n = pCell.phenotype.cycle.data.currentPhaseIndex;

        // this multiplier is for linear interpolation of the oxygen value 
        double multiplier = 1.0;
        if( pO2 < pCell.parameters.o2_proliferation_saturation )
        {
            multiplier = ( pO2 - pCell.parameters.o2_proliferation_threshold )
                    / ( pCell.parameters.o2_proliferation_saturation - pCell.parameters.o2_proliferation_threshold );
        }
        if( pO2 < pCell.parameters.o2_proliferation_threshold )
        {
            multiplier = 0.0;
        }

        // now, update the appropriate cycle transition rate 
        //            pCell.parameters.pReference_live_phenotype.cycle.data.transition_rate( start_phase_index, end_phase_index )\
        //            double base = pCell.parameters.pReference_live_phenotype.cycle.data.transition_rate( start_phase_index, end_phase_index );
        //            System.out.println( pCell.toString() + " " + multiplier );
        phenotype.cycle.data.modifyTransitionRate( start_phase_index, end_phase_index, multiplier );

        // Update necrosis rate
        multiplier = 0.0;
        if( pO2 < pCell.parameters.o2_necrosis_threshold )
        {
            multiplier = ( pCell.parameters.o2_necrosis_threshold - pO2 )
                    / ( pCell.parameters.o2_necrosis_threshold - pCell.parameters.o2_necrosis_max );
        }
        if( pO2 < pCell.parameters.o2_necrosis_max )
        {
            multiplier = 1.0;
        }

        // now, update the necrosis rate 
        pCell.phenotype.death.rates.set( necrosis_index, multiplier * pCell.parameters.max_necrosis_rate );
        // check for deterministic necrosis 
        if( pCell.parameters.necrosis_type == PhysiCellConstants.deterministic_necrosis && multiplier > 1e-16 )
        {
            pCell.phenotype.death.rates.set( necrosis_index, 9e99 );
        }
    }
}