package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.PhaseEntry;
import ru.biosoft.physicell.core.Phenotype;

public class StandardLivePhaseEntry implements PhaseEntry
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        // the cell wants to double its volume 
        phenotype.volume.target_solid_nuclear *= 2.0;
        phenotype.volume.target_solid_cytoplasmic *= 2.0;
    }
}