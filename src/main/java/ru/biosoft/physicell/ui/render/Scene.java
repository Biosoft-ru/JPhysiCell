package ru.biosoft.physicell.ui.render;

import java.util.ArrayList;
import java.util.List;

public class Scene
{
    private List<Mesh> spheres = new ArrayList<Mesh>();
    
    private List<Mesh> circles = new ArrayList<Mesh>();
    
    public void addCircle(Mesh mesh)
    {
        circles.add( mesh );
    }
    
    public void addSphere(Mesh mesh)
    {
        spheres.add( mesh );
    }

    public void clearCircles()
    {
        circles.clear();
    }
    
    public void clear()
    {
        spheres.clear();
        circles.clear();
    }

    public int getSpheresCount()
    {
        return spheres.size();
    }

    public Iterable<Mesh> getSpheres()
    {
        return spheres;
    }

    public Iterable<Mesh> getCircles()
    {
        return circles;
    }

    
    public void sort()
    {
        spheres.sort( null );
    }
}
