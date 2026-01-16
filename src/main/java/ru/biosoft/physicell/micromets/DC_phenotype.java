package ru.biosoft.physicell.micromets;

import java.util.Set;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.standard.StandardModels;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Model;

public class DC_phenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        Model model = pCell.getModel();
        CellDefinition pCD = model.getCellDefinition( "DC" );
        int apoptosis_index = phenotype.death.findDeathModelIndex( "Apoptosis" );

        // no apoptosis until activation (resident DC in constant number for homeostasis)
        if( pCell.customData.get( "activated_immune_cell" ) < 0.5 )
        {
            phenotype.death.rates.set( apoptosis_index, 0.0 );
        }
        else
        {
            phenotype.death.rates.set( apoptosis_index, pCD.phenotype.death.rates.get( apoptosis_index ) );
        }

        // (Adrianne) get type of CD8+ T cell
        int CD8_Tcell_type = model.getCellDefinition( "CD8 Tcell" ).type;
        Cell pTempCell = null;

        if( pCell.customData.get( "activated_immune_cell" ) > 0.5 ) // (Adrianne) activated DCs that don't leave the tissue can further activate CD8s increasing their proliferation rate and attachment rates
        {

            Set<Cell> neighbors = pCell.cells_in_my_container(); // (Adrianne) find cells in a neighbourhood of DCs
            //            int n = 0;
            //            Cell pTestCell = neighbors[n];
            //            while( n < neighbors.size() )
            //            {
            //                pTestCell = neighbors[n];
            for( Cell pTestCell : neighbors )
            {
                // (Adrianne) find the euclidean distance between the DC and the cell it's testing
                double cell_cell_distance = Math
                        .sqrt( ( pTestCell.position[0] - pCell.position[0] ) * ( pTestCell.position[0] - pCell.position[0] )
                                + ( pTestCell.position[1] - pCell.position[1] ) * ( pTestCell.position[1] - pCell.position[1] ) );
                double radius_DC = pCell.phenotype.geometry.radius; // (Adrianne) radius of DC)
                double radius_test_cell = pTestCell.phenotype.geometry.radius; // (Adrianne) radius of test cell)

                // (Adrianne) check if any neighbour cells are live T cells and that they are close enough to the DC
                if( pTestCell != pCell && pTestCell.phenotype.death.dead == false && pTestCell.type == CD8_Tcell_type
                        && cell_cell_distance <= model.getParameterDouble( "epsilon_distance" ) * ( radius_DC + radius_test_cell ) )
                {

                    pTestCell.customData.set( "cell_attachment_rate", model.getParameterDouble( "DC_induced_CD8_attachment" ) ); // (Adrianne) DC induced T cell attachement rate

                    // (Adrianne) finding the G0G1 and S phase index and setting the transition rate to be non zero so that CD8 T cells start proliferating after interacting with DC
                    int cycle_G0G1_index = StandardModels.flow_cytometry_separated_cycle_model
                            .findPhaseIndex( PhysiCellConstants.G0G1_phase );
                    int cycle_S_index = StandardModels.flow_cytometry_separated_cycle_model.findPhaseIndex( PhysiCellConstants.S_phase );
                    pTestCell.phenotype.cycle.data.setTransitionRate( cycle_G0G1_index, cycle_S_index,
                            model.getParameterDouble( "DC_induced_CD8_proliferation" ) );

                    //                    n = neighbors.size();
                    break;
                }

            }
            return;
        }
        else
        {
            //Activation of DC cells - attach with death cancer cells
            // if this returns non-NULL, we're now attached to at least one cell
            if( Util.immune_cell_check_neighbors_for_attachment( pCell, dt ) != null ) // Look around by cancer or lung cells to attach
            {
                //                for (int i = 0; i < pCell.state.numberAttachedCells(); i++){
                //                    pTempCell = pCell.state.attachedCells.get(i);
                //                    #pragma omp critical
                for( Cell cell : pCell.state.attachedCells )
                {
                    pCell.customData.set( "activated_immune_cell", 1.0 );
                    phenotype.motility.isMotile = false;
                }
            }

        }
        return;
    }

}

