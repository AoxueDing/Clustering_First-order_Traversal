package cmsc420_s23;

/**
 * @author Aoxue Ding
 * @Package cmsc420_s23
 * @date 5/3/23 11:37 AM
 */
public class AssignedPair<LPoint extends LabeledPoint2D> {
    LPoint site; // a site
    LPoint center; // its assigned center
    double distanceSq; // the squared distance between them
    public int compareTo(AssignedPair o) {
        if (this.distanceSq < o.distanceSq) {
            return -1;
        } else if (this.distanceSq > o.distanceSq) {
            return 1;
        } else {
            return 0;
        }
    }
}
