package ru.biosoft.physicell.covid;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.standard.StandardModels;

public class CD4_Tcell_phenotype extends UpdatePhenotype
{
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        int cycle_G0G1_index = StandardModels.flow_cytometry_separated_cycle_model.findPhaseIndex( PhysiCellConstants.G0G1_phase );
        int cycle_S_index = StandardModels.flow_cytometry_separated_cycle_model.findPhaseIndex( PhysiCellConstants.S_phase );
        int virus_index = pCell.getMicroenvironment().findDensityIndex( "virion" );
//        int nV_external = virus_index;
//        double virus_amount = pCell.nearest_density_vector()[virus_index];


        int apoptosis_index = pCell.phenotype.death.findDeathModelIndex( "apoptosis" );

        // (AJ-V5) Model for T cell proliferation and death
        int generation_value = (int)pCell.customData.get( "generation" );

        if( pCell.phenotype.cycle.data.elapsedTimePhase < 6 && pCell.phenotype.cycle.data.currentPhaseIndex == 0 )
        {
            pCell.customData.set( "generation", pCell.customData.get( "generation" ) - 1 );
        }
        if( generation_value < 0 )
        {
            pCell.phenotype.death.rates.set( apoptosis_index, pCell.getModel().getParameterDouble( "Death_rates_of_old_Tcells" ) );// [apoptosis_index] = pCell.getModel().getParameterDouble( "Death_rates_of_old_Tcells" );
            pCell.phenotype.death.rates.set( apoptosis_index, 100d ); // new death rate of T cells when they have exceeded generation
            pCell.phenotype.cycle.data.setTransitionRate( cycle_G0G1_index, cycle_S_index, 0 );//  transitionRate ( cycle_G0G1_index, cycle_S_index ) = 0;
        }
    }
}