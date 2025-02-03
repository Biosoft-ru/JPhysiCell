package ru.biosoft.physicell.sample_projects.cancer_metabolism;

import java.awt.Color;

import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.ui.AgentColorer;

public class CancerMetabolismVisualizer implements AgentColorer
{
    @Override
    public Color[] findColors(Cell cell)
    {
        if( !cell.phenotype.death.dead )
            return new Color[] {Color.black, Color.black};

        return new Color[] {Color.blue, Color.black};
    }
}