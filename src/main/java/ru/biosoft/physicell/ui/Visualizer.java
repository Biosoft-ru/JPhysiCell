package ru.biosoft.physicell.ui;

import java.awt.image.BufferedImage;
import java.io.IOException;
import ru.biosoft.physicell.biofvm.Microenvironment;


public interface Visualizer
{
    public String getName();

    public void addResultGenerator(ResultGenerator generator);
    public void clearResultGenerators();
    public Iterable<ResultGenerator> getGenerators();
    public void setAgentColorer(AgentColorer colorer);

    public void init() throws IOException;

    public void saveResult(Microenvironment m, double t) throws IOException;

    public void update(BufferedImage image) throws IOException;
    public void finish() throws IOException;

    public BufferedImage getImage(Microenvironment m, double t) throws Exception;
}