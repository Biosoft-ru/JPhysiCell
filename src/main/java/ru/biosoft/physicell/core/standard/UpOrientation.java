package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.set_orientation;
import ru.biosoft.physicell.core.Phenotype;

public class UpOrientation extends set_orientation
{
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        pCell.state.orientation[0] = 0.0;
        pCell.state.orientation[1] = 0.0;
        pCell.state.orientation[2] = 1.0;
    }
}
