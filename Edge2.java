package algorithms;

import java.awt.Point;


public class Edge implements Comparable<Edge>{
	Point pointA;
    Point pointB;
    int indexA;
    int indexB;
    double dist;

    public Edge(Point pointA, Point pointB, int indexA, int indexB, double dist) {
        this.pointA = pointA;
        this.pointB = pointB;
        this.indexA= indexA;
        this.indexB = indexB;
        this.dist = dist;
    }
    
    public double distance(){ 
    	return dist;
    }

	@Override
	public int compareTo(Edge o) {
		if (this.dist > o.distance())
			return 1;
		else if (this.dist < o.distance())
			return -1;
		else
			return 0;
	}

}
