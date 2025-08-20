use std::ops::{Add, Mul, MulAssign};
use crate::matrix::Matrix;

impl<T> Matrix<T> {
    /// Multiplies every `self[i, j]` with `coefficient_row[0, j]`
    fn seep(&mut self, coefficient_row: Matrix<T>) where T: MulAssign + Clone {
        assert_eq!(coefficient_row.rows, 1);
        assert_eq!(self.columns, coefficient_row.columns);

        for row in 0..self.rows {
            for column in 0..self.columns {
                self[(row, column)] *= coefficient_row[(0, column)].clone();
            }
        }
    }

    /// Turns `self` into a column of the row sums
    fn squish(&mut self) where T: Add<Output = T> {
        assert!(self.columns > 0);

        self.columns = 1;

        for row in &mut self.cells {
            let old_row = std::mem::take(row);
            let row_sum = old_row.into_iter().reduce(T::add).unwrap();
            *row = vec![row_sum];
        }
    }
}

impl<T> Mul for Matrix<T> where T: Add<Output = T> + MulAssign + Clone + Default {
    type Output = Matrix<T>;

    fn mul(self, other: Self) -> Self::Output {
        let a = self;
        let mut b = other;

        b.transpose();

        let mut c = Self::new(a.rows, 0);

        for b_row in 0..b.rows {
            let mut a_clone = a.clone();
            let b_row_matrix = b.row(b_row);
            a_clone.seep(b_row_matrix);
            a_clone.squish();
            c.attach_to_the_right(a_clone);
        }

        c
    }
}

#[cfg(test)]
mod tests {
    use crate::matrix::Matrix;

    #[test]
    fn test_1() {
        let a = Matrix::from_cells([
            [1, 0, 1],
            [0, 1, 2],
        ]);
        let b = Matrix::from_cells([
            [1, 1],
            [1, 2],
            [0, 1],
        ]);
        let c = a * b;
        let expected_c = Matrix::from_cells([
            [1, 2],
            [1, 4],
        ]);
        assert_eq!(c, expected_c);
    }

    #[test]
    fn test_2() {
        let a = Matrix::from_cells([
            [1, 2],
            [4, 5],
        ]);
        let b = Matrix::from_cells([
            [5, 6],
            [8, 9],
        ]);
        let c = a * b;
        let expected_c = Matrix::from_cells([
            [21, 24],
            [60, 69],
        ]);
        assert_eq!(c, expected_c);
    }


    #[test]
    fn test_3() {
        let a = Matrix::from_cells([
            [0.9, 0.1],
            [0.2, 0.8],
        ]);
        let mut b = a.clone();
        for _ in 0..1000 {
            b = b * a.clone();
        }
        /*
        b = [
            [0.4, 0.6],
            [0.4, 0.6],
        ]
         */
        println!("{b:?}");
    }
}
