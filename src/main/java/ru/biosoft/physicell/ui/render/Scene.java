package ru.biosoft.physicell.ui.render;

import java.util.ArrayList;
import java.util.List;

public class Scene
{
    private List<Mesh> meshes = new ArrayList<Mesh>();

    public void add(Mesh mesh)
    {
        meshes.add( mesh );
    }

    public void clear()
    {
        meshes.clear();
    }
    
    public int getMehesCount()
    {
        return meshes.size();
    }
    
    public Iterable<Mesh> getMeshes()
    {
        return meshes;
    }
    
}
