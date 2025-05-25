<?php
$key = 'tajnykluczszyfr!';
$iv = '1234567890abcdef';

$db = new SQLite3('responses.db');

if ($_SERVER['REQUEST_METHOD'] === 'POST' && isset($_POST['clear'])) {
    $db->exec('DELETE FROM responses');
    echo "<p style='color: red;'>Baza danych została wyczyszczona.</p>";
}

$result = $db->query('SELECT name, surname, age, message, time FROM responses');
?>

<h1>Odpowiedzi:</h1>
<form method="POST">
  <input type="submit" name="clear" value="Wyczyść bazę danych" onclick="return confirm('Na pewno?');">
</form>
<br>
<table border="1" cellpadding="5">
  <tr>
    <th>Imię</th>
    <th>Nazwisko</th>
    <th>Wiek</th>
    <th>Wiadomość</th>
    <th>Czas</th>
  </tr>
  <?php while ($row = $result->fetchArray(SQLITE3_ASSOC)): ?>
    <tr>
      <td><?= htmlspecialchars($row['name']) ?></td>
      <td><?= htmlspecialchars($row['surname']) ?></td>
      <td><?= htmlspecialchars($row['age']) ?></td>
      <td><?= htmlspecialchars(openssl_decrypt($row['message'], 'AES-128-CBC', $key, 0, $iv)) ?></td>
      <td><?= htmlspecialchars($row['time']) ?></td>
    </tr>
  <?php endwhile; ?>
</table>
