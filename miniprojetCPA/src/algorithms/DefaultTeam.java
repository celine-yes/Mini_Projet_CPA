package algorithms;

import java.awt.Point;
import java.util.ArrayList;

public class DefaultTeam {
	
	
	//algo floyd-warshall
	  public int[][] calculShortestPaths(ArrayList<Point> points, int edgeThreshold) {
		int n = points.size();
	    int[][] paths=new int[n][n];
	    n = paths.length;
	    for (int i=0;i<paths.length;i++) for (int j=0;j<paths.length;j++) paths[i][j]=i;
	    
	    double[][] dist= new double[n][n];
	    
	 // Initialisation des distances avec les distances directes entre les points
	    for (int i = 0; i < n; i++) {
	        for (int j = 0; j < n; j++) {
	            Point p1 = points.get(i);
	            Point p2 = points.get(j);
	            if (i==j) {dist[i][i]=0; continue;}
	            
	            double distance = p1.distance(p2);
	            if (distance <= edgeThreshold) {
	                dist[i][j] = distance;
	            } else {
	                dist[i][j] = Double.POSITIVE_INFINITY;
	            }
	            paths[i][j]=j;
	        }
	    }
	    
	    for (int k=0; k<n ; k++) {
	    	for(int i=0; i<n ; i++) {
	    		for(int j=0; j<n ; j++) {
	    			if (dist[i][k] + dist[k][j] < dist[i][j]) {
	                    dist[i][j] = dist[i][k] + dist[k][j];
	                    paths[i][j]=paths[i][k];
	                }
	    		}
	    	}
	    }
	    return paths;
	  }
	  
	  
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		  int n = points.size();
		  int m = hitPoints.size();

		  //graphe pondéré complet K
		  ArrayList<Edge> K = new ArrayList<Edge>();
		  for (int i = 0; i < m; i++) {
		      for (int j = 0; j < m; j++) {
		          int indexI = points.indexOf(hitPoints.get(i));
		          int indexJ = points.indexOf(hitPoints.get(j));
		          K.add(new Edge(points.get(indexI), points.get(indexJ), indexI, indexJ)) ;
		      }
		  }
		 
		  //T0 :arbre couvrant à partir de kruskal sur K
		  K = sort(K);

		  ArrayList<Edge> kruskal = new ArrayList<Edge>();
		  Edge current;
		  LabelPoint forest = new LabelPoint(points);
		  while (K.size()!=0) {
		    current = K.remove(0);
		    if (forest.getLabel(current.pointA)!=forest.getLabel(current.pointB)) {
		      kruskal.add(current);
		      forest.setLabel(forest.getLabel(current.pointA),forest.getLabel(current.pointB));
		    }
		  }
		  
		  //Tree2D t0 = edgesToTree(kruskal,kruskal.get(0).pointA);
		  
		  //Construction de H à partir de t0
		  int [][] path = calculShortestPaths(points, edgeThreshold);
		  ArrayList<Edge> H = new ArrayList<Edge>();
		  
		  for(int i=0; i<kruskal.size(); i++) {
			  int indexStart = kruskal.get(i).indexA;
			  int indNextPoint = path[indexStart][kruskal.get(i).indexB];
			  while(indNextPoint != indexStart) {
				  H.add(new Edge(points.get(indexStart), points.get(indNextPoint), indexStart, indNextPoint));
				  indexStart = indNextPoint;
				  indNextPoint = path[indexStart][kruskal.get(i).indexB];
				  
			  }
		  }
		 
	    return edgesToTree(H,H.get(0).pointA);
	  }
	  
	  
	  private ArrayList<Edge> sort(ArrayList<Edge> edges) {
		    if (edges.size()==1) return edges;

		    ArrayList<Edge> left = new ArrayList<Edge>();
		    ArrayList<Edge> right = new ArrayList<Edge>();
		    int n=edges.size();
		    for (int i=0;i<n/2;i++) { left.add(edges.remove(0)); }
		    while (edges.size()!=0) { right.add(edges.remove(0)); }
		    left = sort(left);
		    right = sort(right);

		    ArrayList<Edge> result = new ArrayList<Edge>();
		    while (left.size()!=0 || right.size()!=0) {
		      if (left.size()==0) { result.add(right.remove(0)); continue; }
		      if (right.size()==0) { result.add(left.remove(0)); continue; }
		      if (left.get(0).distance() < right.get(0).distance()) result.add(left.remove(0));
		      else result.add(right.remove(0));
		    }
		    return result;
	  }
	  
	  private Tree2D edgesToTree(ArrayList<Edge> edges, Point root) {
		    ArrayList<Edge> remainder = new ArrayList<Edge>();
		    ArrayList<Point> subTreeRoots = new ArrayList<Point>();
		    Edge current;
		    while (edges.size()!=0) {
		      current = edges.remove(0);
		      if (current.pointA.equals(root)) {
		        subTreeRoots.add(current.pointB);
		      } else {
		        if (current.pointB.equals(root)) {
		          subTreeRoots.add(current.pointA);
		        } else {
		          remainder.add(current);
		        }
		      }
		    }

		    ArrayList<Tree2D> subTrees = new ArrayList<Tree2D>();
		    for (Point subTreeRoot: subTreeRoots) subTrees.add(edgesToTree((ArrayList<Edge>)remainder.clone(),subTreeRoot));

		    return new Tree2D(root, subTrees);
		  }
	  
	  
	  ////////////////////////////////////////////////////////////////////
	  
	  
	  
  public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
	  final double budget = 1664.0; // Budget fixé
	  int n = points.size();
	  int m = hitPoints.size();

	  //graphe pondéré complet K
	  ArrayList<Edge> K = new ArrayList<Edge>();
	  for (int i = 0; i < m; i++) {
	      for (int j = 0; j < m; j++) {
	          int indexI = points.indexOf(hitPoints.get(i));
	          int indexJ = points.indexOf(hitPoints.get(j));
	          K.add(new Edge(points.get(indexI), points.get(indexJ), indexI, indexJ)) ;
	      }
	  }
	 
	  K = sort(K);

	  ArrayList<Edge> kruskal = new ArrayList<Edge>();
	  Edge current;
	  LabelPoint forest = new LabelPoint(points);
	  while (K.size()!=0) {
	    current = K.remove(0);
	    if (forest.getLabel(current.pointA)!=forest.getLabel(current.pointB)) {
	      kruskal.add(current);
	      forest.setLabel(forest.getLabel(current.pointA),forest.getLabel(current.pointB));
	    }
	  }
	  
	 kruskal = sort(kruskal);
	 
	 
	// Calcul du chemin le plus court pour chaque paire de points
	
	    int[][] paths=new int[n][n];
	    n = paths.length;
	    for (int i=0;i<paths.length;i++) for (int j=0;j<paths.length;j++) paths[i][j]=i;
	    
	    double[][] dist= new double[n][n];
	    
	 // Initialisation des distances avec les distances directes entre les points
	for (int i = 0; i < n; i++) {
	    for (int j = 0; j < n; j++) {
	        Point p1 = points.get(i);
	        Point p2 = points.get(j);
	        if (i==j) {dist[i][i]=0; continue;}
	        
	        double distance = p1.distance(p2);
	        if (distance <= edgeThreshold) {
	            dist[i][j] = distance;
	        } else {
	            dist[i][j] = Double.POSITIVE_INFINITY;
	        }
	        paths[i][j]=j;
	    }
	}
	
	for (int k=0; k<n ; k++) {
		for(int i=0; i<n ; i++) {
			for(int j=0; j<n ; j++) {
				if (dist[i][k] + dist[k][j] < dist[i][j]) {
	                dist[i][j] = dist[i][k] + dist[k][j];
	                paths[i][j]=paths[i][k];
	            }
			}
		}
	}
	 
	 
	  ArrayList<Edge> H = new ArrayList<Edge>();
	  double currentBudget = budget; // Initialise le budget restant
	  Point firstHPoint = hitPoints.get(0);
	  Point firstPoint = firstHPoint;
	  ArrayList<Point> remainingHitPoints = new ArrayList<Point>(hitPoints);
	  ArrayList<Point> addedHitPoints = new ArrayList<Point>();
	  remainingHitPoints.remove(firstPoint); // Enlever le premier hitPoint déjà traité
	  addedHitPoints.add(firstPoint);
	  
	  // Continuer tant qu'il reste du budget et qu'il y a des hitPoints à traiter
	  while (currentBudget > 0 && !remainingHitPoints.isEmpty()) {
	      Edge shortestEdge = null;
	      double shortestDistance = Double.POSITIVE_INFINITY;
	      Point endPoint = null;
	      firstPoint = null;
	      for (Point fPoint : addedHitPoints) {
	    	  for (Point hitPoint : remainingHitPoints) {
	    		  int indexStart = points.indexOf(fPoint);
		    	  int indexEnd = points.indexOf(hitPoint);
		    	  double distance = dist[indexStart][indexEnd];
		    	  if(distance != indexStart 
		    			  && distance < shortestDistance
		    			  && distance <= currentBudget) {
		    		  shortestDistance = distance;
	    	          shortestEdge = new Edge(firstPoint, hitPoint, points.indexOf(fPoint), indexEnd);
		    		  endPoint = hitPoint;
		    		  firstPoint = fPoint;
		    	  }
		    	    
	    	  }
	      }

	      // Ajouter le chemin le plus court à H si trouvé
	      if (shortestEdge != null) {
	          currentBudget -= shortestDistance; // Mettre à jour le budget
	          remainingHitPoints.remove(endPoint); // Enlever le hitPoint connecté
	          addedHitPoints.add(endPoint);

	          // Ajouter toutes les arêtes du chemin à H
	          int currentIndex = shortestEdge.indexA;
	          int endIndex = shortestEdge.indexB;
	          while (currentIndex != endIndex) {
	              int nextIndex = paths[currentIndex][endIndex];
	              H.add(new Edge(points.get(currentIndex), points.get(nextIndex), currentIndex, nextIndex));
	              currentIndex = nextIndex;
	          }

	          firstPoint = endPoint; // Mettre à jour le point de départ pour la prochaine itération
	      } else {
	          // Aucun chemin trouvé dans le budget restant, arrêter la boucle
	          break;
	      }
	  }

	  // Convertir H en Tree2D
	  return edgesToTree(H, hitPoints.get(0));
  }
  
}
