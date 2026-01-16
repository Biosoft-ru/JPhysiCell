package ru.biosoft.physicell.micromets;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.standard.StandardModels;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;

public class CD4_Tcell_plenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        int cycle_G0G1_index =  StandardModels.flow_cytometry_separated_cycle_model.findPhaseIndex( PhysiCellConstants.G0G1_phase );
        int cycle_S_index = StandardModels.flow_cytometry_separated_cycle_model.findPhaseIndex( PhysiCellConstants.S_phase );

        int apoptosis_index = pCell.phenotype.death.findDeathModelIndex( "apoptosis" );
        int generation_index = pCell.customData.findVariableIndex( "division_generation" );

        // Model of proliferation based on generation
        if( pCell.phenotype.cycle.data.elapsedTimePhase < 6 && pCell.phenotype.cycle.data.currentPhaseIndex == 0 )
        {
            pCell.customData.set( generation_index, pCell.customData.get( generation_index ) + 1 );
        }

        if( pCell.customData.get(generation_index) > 4 )
        {
            pCell.phenotype.death.rates.set(apoptosis_index, 100.0); // new death rate of T cells when they have exceeded generation
            pCell.phenotype.cycle.data.setTransitionRate(  cycle_G0G1_index, cycle_S_index , 0);
        }

        return;
    }
}
