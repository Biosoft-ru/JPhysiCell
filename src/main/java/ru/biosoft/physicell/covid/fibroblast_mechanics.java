package ru.biosoft.physicell.covid;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;

public class fibroblast_mechanics extends CustomCellRule
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        if( ModelCovid.check_for_out_of_bounds( pCell, 10.0 ) )
            ( (ModelCovid)pCell.getModel() ).cells_to_move_from_edge.add( pCell );
    }
}