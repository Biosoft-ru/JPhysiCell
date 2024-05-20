package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.core.CellFunctions.Contact;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;
import ru.biosoft.physicell.core.CellFunctions.DistanceCalculator;
import ru.biosoft.physicell.core.CellFunctions.MembraneInteractions;
import ru.biosoft.physicell.core.CellFunctions.UpdateMigrationBias;
import ru.biosoft.physicell.core.CellFunctions.UpdatePhenotype;
import ru.biosoft.physicell.core.CellFunctions.UpdateVelocity;
import ru.biosoft.physicell.core.CellFunctions.VolumeUpdate;
import ru.biosoft.physicell.core.CellFunctions.instantiate_cell;
import ru.biosoft.physicell.core.CellFunctions.set_orientation;
import ru.biosoft.physicell.sample_projects.pred_prey_farmer.AvoidBoundariesRule;
import ru.biosoft.physicell.sample_projects.pred_prey_farmer.WrapBoundariesRule;

public class FunctionRegistry
{
    public static VolumeUpdate[] getVolumeFunctions()
    {
        return new VolumeUpdate[] {new StandardVolumeUpdate()};
    }

    public static UpdatePhenotype[] getUpdatePhenotypeFunctions()
    {
        return new UpdatePhenotype[] {new O2based()};
    }

    public static CustomCellRule[] getCustomRules()
    {
        return new CustomCellRule[] {new AvoidBoundariesRule(), new WrapBoundariesRule()};
    }

    public static UpdateMigrationBias[] getUpdateMigrationFunctions()
    {
        return new UpdateMigrationBias[] {};
    }

    public static UpdateVelocity[] getVelocityFunctions()
    {
        return new UpdateVelocity[] {new StandardUpdateVelocity()};
    }

    public static Contact[] getContactFunctions()
    {
        return new Contact[] {new StandardElasticContact()};
    }

    public static instantiate_cell[] getIntsnatiateFunctions()
    {
        return new instantiate_cell[] {};
    }

    public static set_orientation[] getOrientationFunctions()
    {
        return new set_orientation[] {new UpOrientation()};
    }

    public static MembraneInteractions[] getMembraneInteractionFunctions()
    {
        return new MembraneInteractions[] {new DomainEdgeAvoidance()};
    }

    public static DistanceCalculator[] getDistanceCalculatorFunctions()
    {
        return new DistanceCalculator[] {new DomainEdgeDistance()};
    }
}
