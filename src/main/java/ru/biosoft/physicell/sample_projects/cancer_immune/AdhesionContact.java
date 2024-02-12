package ru.biosoft.physicell.sample_projects.cancer_immune;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.contact_function;
import ru.biosoft.physicell.core.Phenotype;

public class AdhesionContact implements contact_function
{
    public void execute(Cell pActingOn, Phenotype pao, Cell pAttachedTo, Phenotype pat, double dt)
    {
        double[] displacement = VectorUtil.newDiff( pAttachedTo.position, pActingOn.position );
        double maxElasticDisplacement = pao.geometry.radius * pao.mechanics.relative_detachment_distance;
        double maxDisplacementSquared = maxElasticDisplacement * maxElasticDisplacement;

        // detach cells if too far apart 
        if( VectorUtil.norm_squared( displacement ) > maxDisplacementSquared )
        {
            Cell.detach_cells( pActingOn, pAttachedTo );
            return;
        }
        VectorUtil.axpy( pActingOn.velocity, pao.mechanics.attachment_elastic_constant, displacement );
    }
}