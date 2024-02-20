package algorithms;

import java.awt.Point;


public class Edge{
	Point pointA;
    Point pointB;
    int indexA;
    int indexB;

    public Edge(Point pointA, Point pointB, int indexA, int indexB) {
        this.pointA = pointA;
        this.pointB = pointB;
        this.indexA= indexA;
        this.indexB = indexB;
    }
    
    public double distance(){ 
    	return pointA.distance(pointB);
    }

}
