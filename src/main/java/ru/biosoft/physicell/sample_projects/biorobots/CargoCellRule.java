package ru.biosoft.physicell.sample_projects.biorobots;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.SignalBehavior;

public class CargoCellRule extends UpdatePhenotype
{
    private double elasticCoefficient;

    public CargoCellRule(Model model)
    {
        elasticCoefficient = model.getParameterDouble( "elastic_coefficient" );
    }

    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        SignalBehavior.setSingleBehavior( pCell, "cell-cell adhesion elastic constant", elasticCoefficient );
    }
}