package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;

public class O2based extends UpdatePhenotype
{
    int start_phase_index = 0; // Q_phase_index; 
    int end_phase_index = 0; // K_phase_index;
    int necrosis_index = 0;
    boolean indicesInitiated = false;

    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        if( phenotype.death.dead == true )
            return;

        int oxygenIndex = pCell.getMicroenvironment().findDensityIndex( "oxygen" );

        if( !indicesInitiated )
            init( phenotype );

        if( !indicesInitiated )// don't continue if we never "figured out" the current cycle model. 
            return;

        // sample the microenvironment to get the pO2 value 
        double pO2 = ( pCell.nearest_density_vector() )[oxygenIndex]; // PhysiCellConstants.oxygen_index]; 

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


    // set up shortcuts to find the Q and K(1) phases (assuming Ki67 basic or advanced model)
    private void init(Phenotype phenotype)
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
                indicesInitiated = true;
            }
            if( phenotype.cycle.code == PhysiCellConstants.advanced_Ki67_cycle_model )
            {
                end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.Ki67_positive_premitotic );
                indicesInitiated = true;
            }
        }

        // live model 
        if( phenotype.cycle.code == PhysiCellConstants.live_cells_cycle_model )
        {
            start_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.live );
            necrosis_index = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );
            end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.live );
            indicesInitiated = true;
        }

        // cytometry models 
        if( phenotype.cycle.code == PhysiCellConstants.flow_cytometry_cycle_model
                || phenotype.cycle.code == PhysiCellConstants.flow_cytometry_separated_cycle_model )
        {
            start_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.G0G1_phase );
            necrosis_index = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );
            end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.S_phase );
            indicesInitiated = true;
        }

        if( phenotype.cycle.code == PhysiCellConstants.cycling_quiescent_model )
        {
            start_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.quiescent );
            necrosis_index = phenotype.death.findDeathModelIndex( PhysiCellConstants.necrosis_death_model );
            end_phase_index = phenotype.cycle.findPhaseIndex( PhysiCellConstants.cycling );
            indicesInitiated = true;
        }
    }

    @Override
    public String display()
    {
        return "Default O2-based phenotype: cell division and necrosis are based on oxygen density";
    }
}