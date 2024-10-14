package ru.biosoft.physicell.ui;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageTypeSpecifier;
import javax.imageio.ImageWriter;
import javax.imageio.metadata.IIOInvalidTreeException;
import javax.imageio.metadata.IIOMetadata;
import javax.imageio.metadata.IIOMetadataNode;
import javax.imageio.stream.ImageOutputStream;

public class GIFGenerator extends ResultGenerator
{
    private ImageWriter writer = null;
    private ImageOutputStream ios = null;
    private IIOMetadata metadata;

    public GIFGenerator(String folder, String name)
    {
        this( new File( folder + "/" + name ) );
    }

    public GIFGenerator(File file)
    {
        super( file );
    }

    public void init() throws IOException
    {
        writer = ImageIO.getImageWritersByFormatName( "GIF" ).next();
        ios = ImageIO.createImageOutputStream( result );

        ImageTypeSpecifier imageTypeSpecifier = ImageTypeSpecifier.createFromBufferedImageType( BufferedImage.TYPE_INT_ARGB );
        metadata = writer.getDefaultImageMetadata( imageTypeSpecifier, writer.getDefaultWriteParam() );
        configureRootMetadata( 100, true );

        writer.setOutput( ios );
        writer.prepareWriteSequence( null );
    }

    public void update(BufferedImage image) throws IOException
    {
        writer.writeToSequence( new IIOImage( image, null, metadata ), writer.getDefaultWriteParam() );
    }

    public void finish() throws IOException
    {
        writer.endWriteSequence();
        writer.reset();
        writer.dispose();
        ios.flush();
        ios.close();
    }

    private void configureRootMetadata(int delay, boolean loop) throws IIOInvalidTreeException
    {
        String metaFormatName = metadata.getNativeMetadataFormatName();
        IIOMetadataNode root = (IIOMetadataNode)metadata.getAsTree( metaFormatName );

        IIOMetadataNode graphicsControlExtensionNode = getNode( root, "GraphicControlExtension" );
        graphicsControlExtensionNode.setAttribute( "disposalMethod", "restoreToBackgroundColor" );
        graphicsControlExtensionNode.setAttribute( "userInputFlag", "FALSE" );
        graphicsControlExtensionNode.setAttribute( "transparentColorFlag", "FALSE" );
        graphicsControlExtensionNode.setAttribute( "delayTime", Integer.toString( delay / 10 ) );
        graphicsControlExtensionNode.setAttribute( "transparentColorIndex", "0" );

        IIOMetadataNode commentsNode = getNode( root, "CommentExtensions" );
        commentsNode.setAttribute( "CommentExtension", "Created by BioUML" );

        IIOMetadataNode appExtensionsNode = getNode( root, "ApplicationExtensions" );
        IIOMetadataNode child = new IIOMetadataNode( "ApplicationExtension" );
        child.setAttribute( "applicationID", "NETSCAPE" );
        child.setAttribute( "authenticationCode", "2.0" );

        int loopContinuously = loop ? 0 : 1;
        child.setUserObject( new byte[] {0x1, (byte) ( loopContinuously & 0xFF ), (byte) ( ( loopContinuously >> 8 ) & 0xFF )} );
        appExtensionsNode.appendChild( child );
        metadata.setFromTree( metaFormatName, root );
    }

    private static IIOMetadataNode getNode(IIOMetadataNode rootNode, String nodeName)
    {
        int nNodes = rootNode.getLength();
        for( int i = 0; i < nNodes; i++ )
        {
            if( rootNode.item( i ).getNodeName().equalsIgnoreCase( nodeName ) )
            {
                return (IIOMetadataNode)rootNode.item( i );
            }
        }
        IIOMetadataNode node = new IIOMetadataNode( nodeName );
        rootNode.appendChild( node );
        return ( node );
    }
}
