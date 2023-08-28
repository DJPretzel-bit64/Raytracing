import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Raytracing {
    Vec3 rayColour(Ray r, HittableList world) {
        HitRecord rec = new HitRecord();
        if(world.hit(r, new Interval(0, Double.POSITIVE_INFINITY), rec)) {
            return Vec3.times(Vec3.plus(rec.normal, new Vec3(1.0, 1.0, 1.0)), 0.5);
        }

        Vec3 unitDirection = Vec3.unitVector(r.direction());
        double a = 0.5*(unitDirection.y() + 1.0);
        return Vec3.plus(Vec3.times(new Vec3(1.0, 1.0, 1.0), 1.0 - a), Vec3.times(new Vec3(0.5, 0.7, 1.0), a));
    }

    int rgb(Vec3 colour) {
        int r = (int)(colour.x() * 255);
        int g = (int)(colour.y() * 255);
        int b = (int)(colour.z() * 255);
        int rgb = r;
        rgb = (rgb << 8) + g;
        rgb = (rgb << 8) + b;
        return rgb;
    }

    public Raytracing() {
        System.out.println("Starting...");

        // Image

        double aspectRatio = 16 / 9.0;
        int width = 400;

        // Calculate the image height, and ensure that it's at least 1
        int height = Math.max((int) (width / aspectRatio), 1);

        // World

        HittableList world = new HittableList();
        world.add(new Sphere(new Vec3(0, 0, -1), 0.5));
        world.add(new Sphere(new Vec3(0, -100.5, -1), 100));

        // Camera

        double focalLength = 1.0;
        double viewportHeight = 2.0;
        double viewportWidth = viewportHeight * (width / (double)(height));
        Vec3 cameraCenter = new Vec3(0, 0, 0);

        // Calculate the vectors across the horizontal and down the vertical viewport edges.
        Vec3 viewportU = new Vec3(viewportWidth, 0, 0);
        Vec3 viewportV = new Vec3(0, -viewportHeight, 0);

        // Calculate the horizontal and vertical delta vectors from pixel to pixel.
        Vec3 pixelDeltaU = Vec3.divide(viewportU, width);
        Vec3 pixelDeltaV = Vec3.divide(viewportV, height);

        // Calculate the location of the upper left pixel.
        Vec3 viewportUpperLeft = Vec3.minus(Vec3.minus(Vec3.minus(cameraCenter, new Vec3(0, 0, focalLength)), Vec3.divide(viewportU, 2)), Vec3.divide(viewportV, 2));
        Vec3 pixel00_loc = Vec3.plus(viewportUpperLeft, Vec3.times(Vec3.plus(pixelDeltaU, pixelDeltaV), 0.5));

        //Render

        BufferedImage render = new BufferedImage(width, height, BufferedImage.TYPE_3BYTE_BGR);
        File output = new File("Render.png");


        for(int j = 0; j < height; j++) {
            for(int i = 0; i < width; i++) {
                Vec3 pixelCenter = Vec3.plus(pixel00_loc, Vec3.plus((Vec3.times(pixelDeltaU, i)), Vec3.times(pixelDeltaV, j)));
                Vec3 rayDirection = Vec3.minus(pixelCenter, cameraCenter);
                Ray r = new Ray(cameraCenter, rayDirection);

                Vec3 pixelColour = rayColour(r, world);
                render.setRGB(i, j, rgb(pixelColour));
            }
        }


        try {
            ImageIO.write(render, "PNG", output);
        } catch(IOException e) {
            System.out.println("Something went very wrong...");
        }

        System.out.println("Done");
    }

    public static void main(String[] args) {
        Raytracing raytracing = new Raytracing();
    }
}

class Vec3 {
    public double[] e;

    public Vec3() {
        this.e = new double[]{ 0, 0, 0 };
    }

    public Vec3(double e0, double e1, double e2) {
        this.e = new double[]{ e0, e1, e2 };
    }

    public double x() {
        return e[0];
    }

    public double y() {
        return e[1];
    }

    public double z() {
        return e[2];
    }

    public Vec3 negate() {
        return new Vec3(-this.e[0], -this.e[1], -this.e[2]);
    }

    public double index(int i) {
        return this.e[i];
    }

    public void plusEquals(Vec3 v) {
        this.e[0] += v.e[0];
        this.e[1] += v.e[1];
        this.e[2] += v.e[2];
    }

    public void timesEquals(double t) {
        this.e[0] *= t;
        this.e[1] *= t;
        this.e[2] *= t;
    }

    public void dividedEquals(double t) {
        this.timesEquals(1 / t);
    }

    public double length() {
        return Math.sqrt(lengthSquared());
    }

    public double lengthSquared() {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

    public static Vec3 plus(Vec3 u, Vec3 v) {
        return new Vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
    }

    public static Vec3 minus(Vec3 u, Vec3 v) {
        return new Vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
    }

    public static Vec3 times(Vec3 v, double t) {
        return new Vec3(v.e[0] * t, v.e[1] * t, v.e[2] * t);
    }

    public static Vec3 divide(Vec3 v, double t) {
        return Vec3.times(v, 1 / t);
    }

    public static double dot(Vec3 u, Vec3 v) {
        return  u.e[0] * v.e[0] +
                u.e[1] * v.e[1] +
                u.e[2] * v.e[2];
    }

    public static Vec3 cross(Vec3 u, Vec3 v) {
        return new Vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                        u.e[2] * v.e[0] - u.e[0] * v.e[2],
                        u.e[0] * v.e[1] - u.e[1] * v.e[0]);
    }

    public static Vec3 unitVector(Vec3 v) {
        return Vec3.divide(v, v.length());
    }
}

record Ray(Vec3 origin, Vec3 direction) {
    public Vec3 at(double t) {
        return Vec3.plus(this.origin, Vec3.times(this.direction, t));
    }
}

class HitRecord {
    public Vec3 p = new Vec3();
    public Vec3 normal = new Vec3();
    public double t = 0;
    public boolean frontFace;

    void setFaceNormal(Ray r, Vec3 outwardNormal) {
        frontFace = Vec3.dot(r.direction(), outwardNormal) < 0;
        normal = frontFace ? outwardNormal : outwardNormal.negate();
    }
}

interface Hittable {
    boolean hit(Ray r, Interval rayT, HitRecord rec);
}

class HittableList implements Hittable{
    List<Hittable> objects = new ArrayList<>();

    HittableList(Hittable object) {
        add(object);
    }

    HittableList(){}

    public void clear() {
        objects.clear();
    }

    public void add(Hittable object) {
        objects.add(object);
    }

    @Override
    public boolean hit(Ray r, Interval rayT, HitRecord rec) {
        HitRecord tempRec = new HitRecord();
        boolean hitAnything = false;
        double closestSoFar = rayT.max;

        for(Hittable object : objects) {
            if(object.hit(r, new Interval(rayT.min, closestSoFar), tempRec)) {
                hitAnything = true;
                closestSoFar = tempRec.t;

                rec.normal = tempRec.normal;
            }
        }

        return hitAnything;
    }
}

class Sphere implements Hittable {
    private final Vec3 center;
    private final double radius;

    public Sphere(Vec3 center, double radius) {
        this.center = center;
        this.radius = radius;
    }

    @Override
    public boolean hit(Ray r, Interval rayT, HitRecord rec) {
        Vec3 oc = Vec3.minus(r.origin(), center);
        double a = r.direction().lengthSquared();
        double halfB = Vec3.dot(oc, r.direction());
        double c = oc.lengthSquared() - radius * radius;
        double discriminant = halfB * halfB - a * c;

        if(discriminant < 0)
            return false;
        double sqrtD = Math.sqrt(discriminant);

        // Find the nearest root that lies in the acceptable range.
        double root = (-halfB - sqrtD) / a;
        if(!rayT.surrounds(root)) {
            root = (-halfB + sqrtD) / a;
            if(!rayT.surrounds(root))
                return false;
        }

        rec.t = root;
        rec.p = r.at(rec.t);
        Vec3 outwardNormal = Vec3.divide(Vec3.minus(rec.p, this.center), this.radius);
        rec.setFaceNormal(r, outwardNormal);

        return true;
    }
}

class Interval {
    double min, max;

    public Interval() {
        min = Double.NEGATIVE_INFINITY;
        max = Double.POSITIVE_INFINITY;
    }

    public Interval(double min, double max) {
        this.min = min;
        this.max = max;
    }

    public boolean contains(double x) {
        return min <= x && x <= max;
    }

    public boolean surrounds(double x) {
        return min < x && x < max;
    }

    static final Interval empty = new Interval(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);
    static final Interval universe = new Interval(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
}