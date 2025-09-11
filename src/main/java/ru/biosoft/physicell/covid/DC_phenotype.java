package ru.biosoft.physicell.covid;

import java.util.Set;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.standard.StandardModels;

public class DC_phenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        int CD8_Tcell_type = pCell.getModel().getCellDefinition( "CD8 Tcell" ).type;

        /* // (Adrianne) if DC is already activated, then check whether it leaves the tissue
        if( pCell.customData.get("activated_immune_cell"] >  0.5 && pCell.getModel().getRNG().UniformRandom()() < 0.002)
        {
        extern double DM; //declare existance of DC lymph
        // (Adrianne) DC leaves the tissue and so we lyse that DC
        std::cout<<"DC leaves tissue"<<std::endl;
        pCell.lyse_cell(); 
        #pragma omp critical 
        { DM++; } // add one
        return;
        
        } */
        if( pCell.customData.get( "activated_immune_cell" ) > 0.5 ) // (Adrianne) activated DCs that don't leave the tissue can further activate CD8s increasing their proliferation rate and attachment rates
        {

            Set<Cell> neighbors = pCell.cells_in_my_container(); // (Adrianne) find cells in a neighbourhood of DCs
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
                        && cell_cell_distance <= pCell.getModel().getParameterDouble( "epsilon_distance" )
                                * ( radius_DC + radius_test_cell ) )
                {

                    pTestCell.customData.set( "cell_attachment_rate", pCell.getModel().getParameterDouble( "DC_induced_CD8_attachment" ) ); // (Adrianne) DC induced T cell attachement rate

                    // (Adrianne) finding the G0G1 and S phase index and setting the transition rate to be non zero so that CD8 T cells start proliferating after interacting with DC
                    int cycle_G0G1_index = StandardModels.flow_cytometry_separated_cycle_model
                            .findPhaseIndex( PhysiCellConstants.G0G1_phase );
                    int cycle_S_index = StandardModels.flow_cytometry_separated_cycle_model.findPhaseIndex( PhysiCellConstants.S_phase );
                    pTestCell.phenotype.cycle.data.setTransitionRate( cycle_G0G1_index, cycle_S_index,
                            pCell.getModel().getParameterDouble( "DC_induced_CD8_proliferation" ) );
                    break;
                }
            }
            return;
        }
        else
        {
            // (adrianne) DCs become activated if there is an infected cell in their neighbour with greater 1 viral protein or if the local amount of virus is greater than 10
            int virus_index = pCell.getMicroenvironment().findDensityIndex( "virion" );
            double virus_amount = pCell.nearest_density_vector()[virus_index];
            if( virus_amount * pCell.getMicroenvironment().mesh.voxels[1].volume > pCell.getModel()
                    .getParameterDouble( "virions_needed_for_DC_activation" ) ) // (Adrianne) amount of virus in local voxel with DC is greater than 10
            {

                pCell.customData.set( "activated_immune_cell", 1.0 ); // (Adrianne) DC becomes activated
                //System.out.println( "Activated virions" );
            }
            else
            {
                for( Cell pTestCell : pCell.cells_in_my_container() )
                {
                    //                    pTestCell = neighbors[n]; 
                    // if it is not me and the target is dead 
                    int nP = pTestCell.customData.findVariableIndex( "viral_protein" );
                    if( pTestCell != pCell && pTestCell.phenotype.death.dead == false && pTestCell.customData.get( nP ) > 1 )
                    {
                        pCell.customData.set( "activated_immune_cell", 1.0 );
                        //                        System.out.println( "Activated protein" );
                        break;
                    }
                }
                return;
            }
        }
    }
}