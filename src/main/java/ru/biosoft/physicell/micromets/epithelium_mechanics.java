package ru.biosoft.physicell.micromets;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;

public class epithelium_mechanics extends CustomCellRule
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );

        // if I'm dead, don't bother
        if( phenotype.death.dead == true )
        {
            // the cell death functions don't automatically turn off custom functions,
            // since those are part of mechanics.
            // remove_all_adhesions( pCell );

            // Let's just fully disable now.
            pCell.functions.customCellRule = null;
            pCell.functions.contact = null;

            phenotype.secretion.secretionRates[debris_index] = pCell.customData.get( "debris_secretion_rate" );
            return;
        }

        // Plastoelastic mechanics
        int spring_constant_index = pCell.customData.findVariableIndex( "spring_constant" );
        int relaxation_constant_index = pCell.customData.findVariableIndex( "mechanical_relaxation_rate" );
        int ECM_attachment_point_index = pCell.customData.findVariableIndex( "ECM_attachment_point" );
        int mechanical_strain_displacement_index = pCell.customData.findVectorVariableIndex( "mechanical_strain_displacement" );

        // first, update the cell's velocity based upon the elastic model
        VectorUtil.axpy( pCell.velocity, pCell.customData.get( spring_constant_index ),
                pCell.customData.vectorVariables.get( mechanical_strain_displacement_index ).value );

//        if( pCell.velocity[0] != 0 || pCell.velocity[1] != 0 )
//            System.out.println( pCell.velocity[0] + " " + pCell.velocity[1] );

        // now, plastic mechanical relaxation
        double plastic_temp_constant = -dt * pCell.customData.get( relaxation_constant_index );
        VectorUtil.axpy( pCell.customData.vectorVariables.get( ECM_attachment_point_index ).value, plastic_temp_constant,
                pCell.customData.vectorVariables.get( mechanical_strain_displacement_index ).value );

        //        axpy( &( pCell->velocity ) , pCell->custom_data[spring_constant_index] , pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value );
        //
        //        // now, plastic mechanical relaxation
        //        static double plastic_temp_constant = -dt * pCell->custom_data[relaxation_constant_index];
        //        axpy( &(pCell->custom_data.vector_variables[ECM_attachment_point_index].value) , plastic_temp_constant , pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value );
        return;
    }

}
