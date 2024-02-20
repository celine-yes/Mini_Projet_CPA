package algorithms;

import java.util.ArrayList;
import java.util.HashMap;
import java.awt.Point;

public class LabelPoint {

	private HashMap<Point, Integer> pointTags;

    public LabelPoint(ArrayList<Point> points) {
        pointTags = new HashMap<>();
        for (int i = 0; i < points.size(); i++) {
            pointTags.put(points.get(i), i);
        }
    }

    public void setLabel(int oldLabel, int newLabel) {
        for (Point p : pointTags.keySet()) {
            if (pointTags.get(p).equals(oldLabel)) {
                pointTags.put(p, newLabel);
            }
        }
    }

    public int getLabel(Point p) {
        return pointTags.getOrDefault(p, 0xBADC0DE);
    }
	
}
