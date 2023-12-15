package ru.biosoft.biofvm.cell;

public abstract class PhaseArrest
{
    public abstract boolean isArrested(Cell pCell, Phenotype phenotype, double dt);

    public static class standard_necrosis_arrest_function extends PhaseArrest
    {
        @Override
        public boolean isArrested(Cell pCell, Phenotype phenotype, double dt)
        {
            // remain in the non-lysed state / phase if volume has not exceeded the 
            // rupture volume 
            if( phenotype.volume.total < phenotype.volume.rupture_volume )
            {
                return true;
            }
            return false;
        }
    }
}