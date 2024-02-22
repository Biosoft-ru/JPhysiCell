package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.CellFunctions.Contact;
import ru.biosoft.physicell.core.Phenotype;

public class StandardElasticContact extends Contact
{
    @Override
    public void execute(Cell pC1, Phenotype p1, Cell pC2, Phenotype p2, double dt)
    {
        if( pC1.position.length != 3 || pC2.position.length != 3 )
        {
            return;
        }

        double[] displacement = VectorUtil.newDiff( pC2.position, pC1.position );
        //        std::vector<double> displacement = pC2.position;
        //        displacement -= pC1.position; 

        // update May 2022 - effective adhesion 
        int ii = CellDefinition.getCellDefinitionIndex( pC1.type );
        int jj = CellDefinition.getCellDefinitionIndex( pC2.type );

        double adhesion_ii = pC1.phenotype.mechanics.attachmentElasticConstant * pC1.phenotype.mechanics.cellAdhesionAffinities[jj];
        double adhesion_jj = pC2.phenotype.mechanics.attachmentElasticConstant * pC2.phenotype.mechanics.cellAdhesionAffinities[ii];

        double effective_attachment_elastic_constant = Math.sqrt( adhesion_ii * adhesion_jj );

        // axpy( &(pC1.velocity) , p1.mechanics.attachment_elastic_constant , displacement ); 
        VectorUtil.axpy( pC1.velocity, effective_attachment_elastic_constant, displacement );
    }
}