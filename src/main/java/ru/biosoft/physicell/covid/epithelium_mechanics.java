package ru.biosoft.physicell.covid;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;

public class epithelium_mechanics extends CustomCellRule
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        //    void epithelium_mechanics(Cell pCell, Phenotype phenotype, double dt)
        //    {
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );

        pCell.isMovable = false;

        // if I'm dead, don't bother 
        if( phenotype.death.dead )
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

        // this is now part of contact_function 
        /*
        // if I'm adhered to something ... 
        if( pCell.state.neighbors.size() > 0 )
        {
        // add the elastic forces 
        extra_elastic_attachment_mechanics( pCell, phenotype, dt );
        }
        */
    }
    
    
}