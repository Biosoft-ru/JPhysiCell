package ru.biosoft.physicell.xml;

import java.io.File;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import ru.biosoft.physicell.biofvm.Microenvironment;

public class ModelReader extends Constants
{

    public void read(File f) throws Exception
    {
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();
        Document doc = builder.parse( f );
        NodeList nodes = doc.getChildNodes();

        Element physicell = findElement( nodes, PHYSICELL_ELEMENT );
        if( physicell == null )
            throw new Exception( "Physicell base element not found" );

        Microenvironment m = new Microenvironment();
        readDomain( physicell, m );
    }

    private void readDomain(Element physicell, Microenvironment m)
    {
        Element domainElement = findElement( physicell, DOMAIN_ELEMENT );
        double minX = 0;
        double maxX = 100;
        double minY = 0;
        double maxY = 100;
        double minZ = 0;
        double maxZ = 100;
        double dx = 10;
        double dy = 10;
        double dz = 10;
        boolean use2D = false;
        NodeList list = domainElement.getChildNodes();
        for( int i = 0; i < list.getLength(); i++ )
        {
            Node node = list.item( i );
            if( node instanceof Element )
            {
                Element el = (Element)node;
                String val = el.getNodeValue();
                switch( ( (Element)node ).getTagName() )
                {
                    case X_MIN:
                        minX = Double.parseDouble( val );
                        break;
                    case X_MAX:
                        maxX = Double.parseDouble( val );
                        break;
                    case Y_MIN:
                        minY = Double.parseDouble( val );
                        break;
                    case Y_MAX:
                        maxY = Double.parseDouble( val );
                        break;
                    case Z_MIN:
                        minZ = Double.parseDouble( val );
                        break;
                    case Z_MAX:
                        maxZ = Double.parseDouble( val );
                        break;
                    case DX:
                        dx = Double.parseDouble( val );
                        break;
                    case DY:
                        dy = Double.parseDouble( val );
                        break;
                    case DZ:
                        dz = Double.parseDouble( val );
                        break;
                    case USE_2D:
                        use2D = Boolean.parseBoolean( val );
                        break;
                }
            }
        }

        m.resizeSpace( minX, maxX, minY, maxY, minZ, maxZ, dx, dy, dz );
    }

    public void readOverall(Element physicell, Microenvironment m)
    {
        double maxTime = 100;
        String maxTimeUnits = "min";
        String timeUnits = "min";
        String spaceUnits = "micron";
        String diffusionUnits = "min";
        String mechanicsUnits = "min";
        String phenotypeUnits = "min";
        double diffusionStep = 0.01;
        double mechanicsStep = 0.1;
        double phenotypeSteps = 6;
        Element overallElement = findElement( physicell, OVERALL );
        NodeList list = overallElement.getChildNodes();
        for( int i = 0; i < list.getLength(); i++ )
        {
            Node node = list.item( i );
            if( node instanceof Element )
            {
                Element el = (Element)node;
                String val = el.getNodeValue();
                switch( ( (Element)node ).getTagName() )
                {
                    case MAX_TIME:
                        maxTime = Double.parseDouble( val );
                        maxTimeUnits = el.getAttribute( UNITS );
                        break;
                    case TIME_UNITS:
                        timeUnits = val;
                        break;
                    case SPACE_UNITS:
                        spaceUnits = val;
                        break;
                    case DT_DIFFUSION:
                        mechanicsUnits = el.getAttribute( UNITS );
                        mechanicsStep = Double.parseDouble( val );
                        break;
                    case DT_MECHANICS:
                        diffusionUnits = el.getAttribute( UNITS );
                        diffusionStep = Double.parseDouble( val );
                        break;
                    case DT_PHENOTYPE:
                        phenotypeUnits = el.getAttribute( UNITS );
                        phenotypeSteps = Double.parseDouble( val );
                        break;

                }
            }
        }
        m.timeUnits = timeUnits;
        m.spatialUnits = spaceUnits;
    }


    @FunctionalInterface
    private static interface ElementParser
    {
        public void parse(Element el, Object obj);
    }

    private Element findElement(Element parent, String tag)
    {
        return findElement( parent.getChildNodes(), tag );
    }

    private Element findElement(NodeList list, String tag)
    {
        for( int i = 0; i < list.getLength(); i++ )
        {
            Node node = list.item( i );
            if( node instanceof Element && ( (Element)node ).getTagName().equals( tag ) )
            {
                return (Element)node;
            }
        }
        return null;
    }

    private void parseElement(Element element, Object obj, ElementParser parser)
    {
        NodeList list = element.getChildNodes();
        for( int i = 0; i < list.getLength(); i++ )
        {
            Node node = list.item( i );
            if( node instanceof Element )
            {
                parser.parse( (Element)node, obj );
            }
        }
    }
}
