package ru.biosoft.physicell.sample_projects.biorobots;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.update_phenotype;

public class DirectorCellRule extends update_phenotype
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        return;
    }
}