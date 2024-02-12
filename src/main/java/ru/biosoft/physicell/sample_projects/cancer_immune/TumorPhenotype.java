package ru.biosoft.physicell.sample_projects.cancer_immune;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.update_phenotype;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.StandardModels;

/**
 * Custom cell phenotype function to scale immunostimulatory factor with hypoxia 
 */
public class TumorPhenotype extends update_phenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        Microenvironment m = pCell.getMicroenvironment();
        int cycle_start_index = StandardModels.live.findPhaseIndex( PhysiCellConstants.live );
        int cycle_end_index = StandardModels.live.findPhaseIndex( PhysiCellConstants.live );
        int oncoprotein_i = pCell.custom_data.find_variable_index( "oncoprotein" );

        // update secretion rates based on hypoxia 
        int o2_index = m.findDensityIndex( "oxygen" );
        int immune_factor_index = m.findDensityIndex( "immunostimulatory factor" );
        double o2 = pCell.nearest_density_vector()[o2_index];

        phenotype.secretion.secretionRates[immune_factor_index] = 10.0;

        new StandardModels.update_cell_and_death_parameters_O2_based().execute( pCell, phenotype, dt );//TODO: use inheritance instead

        // if cell is dead, don't bother with future phenotype changes, set it to secrete the immunostimulatory factor 
        if( phenotype.death.dead == true )
        {
            phenotype.secretion.secretionRates[immune_factor_index] = 10;
            pCell.functions.updatePhenotype = null;
            return;
        }

        // multiply proliferation rate by the oncoprotein 
        //            System.out.println( pCell.toString() + " " + phenotype.cycle.transition_rate( 0, 0 ) );
        double rate = phenotype.cycle.data.getTransitionRate( cycle_start_index, cycle_end_index );
        double factor = pCell.custom_data.get( oncoprotein_i );
        //        System.out.println( pCell.toString() + " " + rate + " " + factor + " " + rate * factor );
        phenotype.cycle.data.setTransitionRate( cycle_start_index, cycle_end_index, rate * factor );
        //            phenotype.cycle.data.modifyTransitionRate( cycle_start_index, cycle_end_index, pCell.custom_data.get( oncoprotein_i ) );
    }
}