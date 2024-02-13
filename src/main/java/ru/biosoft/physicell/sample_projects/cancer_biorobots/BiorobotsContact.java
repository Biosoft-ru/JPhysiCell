package ru.biosoft.physicell.sample_projects.cancer_biorobots;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.contact_function;

public class BiorobotsContact implements contact_function
{
    public void execute(Cell pActingOn, Phenotype pao, Cell pAttachedTo, Phenotype pat, double dt)
    {
        double[] displacement = VectorUtil.newDiff( pAttachedTo.position, pActingOn.position );
        double max_elastic_displacement = pao.geometry.radius * pao.mechanics.relDetachmentDistance;
        double max_displacement_squared = max_elastic_displacement * max_elastic_displacement;
        // detach cells if too far apart
        if( VectorUtil.norm_squared( displacement ) > max_displacement_squared )
        {
            Cell.detachCells( pActingOn, pAttachedTo );
            return;
        }
        VectorUtil.axpy( ( pActingOn.velocity ), pao.mechanics.attachmentElasticConstant, displacement );
    }
}