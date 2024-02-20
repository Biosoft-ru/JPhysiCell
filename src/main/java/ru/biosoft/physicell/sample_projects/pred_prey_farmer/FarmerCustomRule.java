package ru.biosoft.physicell.sample_projects.pred_prey_farmer;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;

/* farmer functions */
public class FarmerCustomRule extends PreyCustomRule
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        super.execute( pCell, phenotype, dt );
    }
}