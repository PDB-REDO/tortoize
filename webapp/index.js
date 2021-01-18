import "core-js/stable";
import "regenerator-runtime/runtime";
import bsCustomFileInput from 'bs-custom-file-input';

class Tortoize {
	constructor(form) {
		this.form = document.forms['tortoize-form'];
		this.table = document.getElementById('tortoize-table');
		this.alert = document.getElementById('tortoize-alert');
		this.form.addEventListener('submit', (evt) => this.handleSubmit(evt));

		this.alert.classList.add('invisible');
	}

	handleSubmit(event) {
		if (event)
			event.preventDefault();
		
		const data = new FormData(this.form);
		this.table.style.display = 'none';
		this.alert.classList.add('invisible');

		let wasOK;
		fetch('tortoize', {
			credentials: 'include',
			method: 'post',
			body: data
		}).then(r => {
			wasOK = r.ok;
			return r.json();
		}).then(r => {
			// console.log(r);
			if (r.model)
				this.process(r.model);
			else if (r.error)
				throw r.error;
			else
				throw 'Reply does not contain data';
		}).catch(err => {
			console.log(err);
			
			this.alert.textContent = `Could not calculate rama z-score: ${err}`;
			this.alert.classList.remove('invisible');
		});
	}

	process(model) {
		const tbody = this.table.querySelector('tbody');
		[...tbody.querySelectorAll('tr')].forEach(row => tbody.removeChild(row));

		for (const [id, data] of Object.entries(model)) {
			const row = document.createElement("tr");
			const tdid = document.createElement("td");
			tdid.textContent = id;
			row.appendChild(tdid);

			for (const f of ['ramachandran-z', 'ramachandran-jackknife-sd', 'torsion-z', 'torsion-jackknife-sd'])
			{
				const td = document.createElement('td');
				td.textContent = (+data[f]).toFixed(2);
				row.appendChild(td);
			}

			tbody.appendChild(row);
		}

		this.table.style.display = 'unset';
	}
}


window.addEventListener('load', () => {
	bsCustomFileInput.init();

	new Tortoize();
});