package ru.biosoft.physicell.micromets;

import java.util.Set;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellConstants;
import ru.biosoft.physicell.core.standard.StandardModels;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Model;

public class epithelium_phenotype extends UpdatePhenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );
        int apoptosis_index = phenotype.death.findDeathModelIndex( "Apoptosis" );

        phenotype.motility.isMotile = false;

        // Mechanical strain
        int strain_index = pCell.customData.findVariableIndex( "mechanical_strain" );
        int ECM_attachment_point_index = pCell.customData.findVectorVariableIndex( "ECM_attachment_point" );
        int mechanical_strain_displacement_index = pCell.customData.findVectorVariableIndex( "mechanical_strain_displacement" );

//        pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value = pCell->custom_data.vector_variables[ECM_attachment_point_index].value;
//        pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value -= pCell->position;
//        pCell->custom_data[strain_index] = norm( pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value );
        
        // Update the strain displacement (ECM - position) - mechanical strain = |ECM - position|
        pCell.customData.vectorVariables.get( mechanical_strain_displacement_index ).value = pCell.customData.vectorVariables
                .get( ECM_attachment_point_index ).value.clone();
        VectorUtil.diff( pCell.customData.vectorVariables.get( mechanical_strain_displacement_index ).value, pCell.position );
        pCell.customData.set( strain_index,
                VectorUtil.norm( pCell.customData.vectorVariables.get( mechanical_strain_displacement_index ).value ) );

//        double[] str =  pCell.customData.vectorVariables.get( mechanical_strain_displacement_index ).value;           
        
        // if I am dead, remove all adhesions
        if( phenotype.death.dead == true )
        {
            // detach all attached cells
            // remove_all_adhesions( pCell );
            phenotype.secretion.secretionRates[debris_index] = pCell.customData.get( "debris_secretion_rate" );
        }

        int cycle_G0G1_index = StandardModels.flow_cytometry_separated_cycle_model.findPhaseIndex( PhysiCellConstants.G0G1_phase );
        int cycle_S_index = StandardModels.flow_cytometry_separated_cycle_model.findPhaseIndex( PhysiCellConstants.S_phase );
        // turn off proliferation
        pCell.phenotype.cycle.data.setTransitionRate( cycle_G0G1_index, cycle_S_index, 0.0 );

        if( strain_based_apoptosis( pCell ) )
        {
            pCell.phenotype.death.rates.set( apoptosis_index, 9e9 );
        }

        // if I am dead, don't bother executing this function again
        if( phenotype.death.dead == true )
        {
            pCell.functions.updatePhenotype = null;
        }

        return;
    }

    boolean strain_based_apoptosis(Cell pCell)
    {
        Model model = pCell.getModel();
        int strain_index = pCell.customData.findVariableIndex( "mechanical_strain" );
        double max_strain = 0.75; // maximum tolerated deformation of lung cell (death) [0.75 microns]
        if( pCell.customData.get( strain_index ) <= max_strain )
            return false;

        Set<Cell> neighbors = pCell.cells_in_my_container();//find cells in neighbourhood
        int cancer_cell_type = pCell.getModel().getCellDefinition( "cancer cell" ).type;
        //        int n = 0;
        //        Cell pTestCell;
        //        while( n < neighbors.size() )
        //        {
        //            pTestCell = neighbors[n];
        // Cell is not himself, live cancer cell
        for( Cell pTestCell : neighbors )
        {
            if( pTestCell != pCell && pTestCell.phenotype.death.dead == false && pTestCell.type == cancer_cell_type )
            {
                double cell_cell_distance = Math
                        .sqrt( ( pTestCell.position[0] - pCell.position[0] ) * ( pTestCell.position[0] - pCell.position[0] )
                                + ( pTestCell.position[1] - pCell.position[1] ) * ( pTestCell.position[1] - pCell.position[1] ) );
                double radius_lung_cell = pCell.phenotype.geometry.radius;
                double radius_cancer_cell = pTestCell.phenotype.geometry.radius;
                if( cell_cell_distance <= model.getParameterDouble( "epsilon_distance" ) * ( radius_lung_cell + radius_cancer_cell ) )
                {
                    return true;
                }
            }
            //            n++;
        }
        return false;
    }
}
