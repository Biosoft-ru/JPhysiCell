package ru.biosoft.physicell.ui;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

public abstract class ResultGenerator
{
    protected File result;

    public ResultGenerator(String folder, String name)
    {
        this( new File( folder + "/" + name ) );
    }

    public ResultGenerator(File file)
    {
        setFile( file );
    }

    public void setFile(File file)
    {
        this.result = file;
    }

    public abstract void init() throws IOException;
    public abstract void update(BufferedImage image) throws IOException;
    public abstract void finish() throws IOException;

    public File getResult()
    {
        return result;
    }
}
