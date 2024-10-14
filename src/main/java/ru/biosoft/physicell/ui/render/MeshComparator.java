package ru.biosoft.physicell.ui.render;

import java.util.Comparator;

public class MeshComparator implements Comparator<Mesh>
{

    @Override
    public int compare(Mesh o1, Mesh o2)
    {
       return Double.compare( o2.center.z, o1.center.z );
    }

}
