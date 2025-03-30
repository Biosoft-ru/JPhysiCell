package ru.biosoft.physicell.ui;

public class ModelData
{
    private Dimension xDim = new Dimension( 0, 100, 10 );
    private Dimension yDim = new Dimension( 0, 100, 10 );
    private Dimension zDim = new Dimension( 0, 100, 10 );
    private boolean use2D = false;
    private String[] substrates;

    public boolean isUse2D()
    {
        return use2D;
    }

    public void setUse2D(boolean use2d)
    {
        use2D = use2d;
    }

    public String[] getSubstrates()
    {
        return substrates;
    }

    public void setSubstrates(String[] substrates)
    {
        this.substrates = substrates;
    }

    public void setXDim(double from, double to, double step)
    {
        xDim = new Dimension( from, to, step );
    }

    public void setYDim(double from, double to, double step)
    {
        yDim = new Dimension( from, to, step );
    }

    public void setZDim(double from, double to, double step)
    {
        zDim = new Dimension( from, to, step );
    }

    public Dimension getXDim()
    {
        return xDim;
    }

    public Dimension getYDim()
    {
        return yDim;
    }

    public Dimension getZDim()
    {
        return zDim;
    }

    public static class Dimension
    {
        private double from;
        private double to;
        private double step;

        public Dimension(double from, double to, double step)
        {
            this.from = from;
            this.to = to;
            this.step = step;
        }

        public double getFrom()
        {
            return from;
        }

        public double getTo()
        {
            return to;
        }

        public double getStep()
        {
            return step;
        }

        public double getLength()
        {
            return to - from;
        }
    }
}
